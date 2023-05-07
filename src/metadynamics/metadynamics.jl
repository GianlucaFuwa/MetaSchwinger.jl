module Metadynamics
	using Base.Threads: nthreads, threadid
	using DelimitedFiles
	using ForwardDiff
	using Polyester
	using StatsBase
	using Optimization 
	using OptimizationOptimJL
	using Zygote

	import ..Gaugefields: Gaugefield
	import ..System_parameters: Params
	import ..Verbose_print: Verbose, Verbose_, println_verbose
	
	struct BiasPotential{O}
		is_static::Bool
		symmetric::Bool

		CVlims::NTuple{2, Float64}
		bin_width::Float64
		w::Float64
		k::Float64

		parametric::Bool
		current_parameters::Union{Nothing, Vector{Float64}}
		lower_bounds::Union{Nothing, Vector{Float64}}
		upper_bounds::Union{Nothing, Vector{Float64}}
		batchsize::Union{Nothing, Int64} 
		cv_storage::Union{Nothing, Vector{Float64}}
		bias_storage::Union{Nothing, Vector{Float64}}
		fullbias_storage::Union{Nothing, Vector{Vector{Float64}}}
		testfun::Union{Nothing, Function}
		minimizer::O

		well_tempered::Bool
		ΔT::Union{Nothing, Float64}

		values::Vector{Float64}
		cv_vals::Vector{Float64}

		exceeded_count::Int64
		fp::Union{Nothing, IOStream}
		KS_fp::Union{Nothing, Verbose}
		
		function BiasPotential(p::Params; instance::Int64 = 1)
			is_static = instance==0 ? true : p.is_static[instance]
			symmetric = p.symmetric

			CVlims = p.CVlims
			bin_width = p.bin_width
			w = p.w
			k = p.k

			parametric = p.parametric

			if parametric == true && is_static == false
				current_parameters = p.potential_parameters
				lower_bounds = p.lower_bounds
				upper_bounds = p.upper_bounds
				batchsize = p.batchsize
				cv_storage = Float64[]
				bias_storage = Float64[]
				fullbias_storage = Vector{Vector{Float64}}(undef, p.Nsweeps)
				testfun = getfield(Metadynamics, Symbol(p.testfun))
				minimizer = getfield(OptimizationOptimJL, Symbol(p.minimizer))()
			else
				current_parameters = nothing
				lower_bounds = nothing
				upper_bounds = nothing
				batchsize = nothing
				cv_storage = nothing
				bias_storage = nothing
				fullbias_storage = nothing
				testfun = nothing
				minimizer = nothing
			end

			well_tempered = p.well_tempered
			ΔT = p.ΔT

			if instance == 0
				values = potential_from_file(p, nothing)
			else
				values = potential_from_file(p, p.usebiases[instance])
			end

			cv_vals = range(CVlims[1], CVlims[2], step=bin_width)

			exceeded_count = 0
			fp = instance==0 ? nothing : open(p.biasfiles[instance], "w")

			if p.parametric == true
				KS_fp = Verbose_(open(
					pwd() * p.measure_dir * "/" * p.testfun * "_" * p.minimizer*".txt",
					"w",
				)) 
			else
				KS_fp = nothing
			end

			return new{typeof(minimizer)}(
			is_static, symmetric,
			CVlims, bin_width, w, k,
			parametric, current_parameters, lower_bounds, upper_bounds, batchsize,
			cv_storage, bias_storage, fullbias_storage, testfun, minimizer,
			well_tempered, ΔT,
			values, cv_vals,
			exceeded_count, fp, KS_fp,
			)
		end	
	end

	function potential_from_file(p::Params, usebias::Union{Nothing,String})
		if usebias === nothing
			len = 1 + round(
				Int,
				(p.CVlims[2] - p.CVlims[1]) / p.bin_width,
				RoundNearestTiesAway
			)
			return zeros(len)
		else
			values = readdlm(usebias, Float64, comments = true)
			len = 1 + round(
				Int,
				(p.CVlims[2] - p.CVlims[1]) / p.bin_width,
				RoundNearestTiesAway
			)
			@assert length(values[:,2]) == len  "Potential length doesn't match parameters"
			return values[:,2]
		end
	end

	function write_to_file!(b::BiasPotential)

		for (idx, cv) in enumerate(b.cv_vals)
			value = b(cv)
			println(b.fp, "$cv $value # Metapotential")
		end

		return nothing
	end

	function Base.flush(b::BiasPotential)
        if b.fp !== nothing
            flush(b.fp)
        end
    end

	function Base.seekstart(b::BiasPotential)
		if b.fp !== nothing
            seekstart(b.fp)
        end
    end

	function Base.setindex!(b::BiasPotential, v, i)
		b.values[min(length(b.values), max(i, 1))] = v
		return nothing
	end

	@inline function Base.getindex(b::BiasPotential,i)
		return b.values[min(length(b.values), max(i, 1))]
	end

	@inline function index(b::BiasPotential, cv)
		grid_index = (cv - b.CVlims[1]) / b.bin_width + 0.5
		return round(Int, grid_index, RoundNearestTiesAway)
	end

	function Base.sum(b::BiasPotential)
		return sum(b.values)
	end

	is_static(b::BiasPotential) = b.is_static
	is_symmetric(b::BiasPotential) = b.symmetric
	is_well_tempered(b::BiasPotential) = b.well_tempered
	is_parametric(b::BiasPotential) = b.parametric

	get_CVlims(b::BiasPotential) = b.CVlims
	get_parameters(b::BiasPotential) = b.current_parameters
	get_bounds(b::BiasPotential) = b.lower_bounds, b.upper_bounds
	get_minimizer(b::BiasPotential) = b.minimizer
	get_batchsize(b::BiasPotential) = b.batchsize
	get_cvstorage(b::BiasPotential) = b.cv_storage
	get_biasstorage(b::BiasPotential) = b.bias_storage
	get_fullbiasstorage(b::BiasPotential) = b.fullbias_storage

	function update_bias!(b::BiasPotential, cv; itrj = nothing)
		if is_static(b)
			return nothing
		elseif is_parametric(b) == false
			update_bias_regular!(b, cv)

			if is_symmetric(b)
				update_bias_regular!(b, -cv)
			end

		elseif is_parametric(b) == true
			update_bias_parametric!(b, cv, itrj)
		end
	end 

	function update_bias_regular!(b::BiasPotential, cv)
		grid_index = index(b, cv)

		if 1 <= grid_index < length(b.values)
			for (idx, current_bin) in enumerate(b.cv_vals)
				fac = b.well_tempered ? exp(-b[idx] / b.ΔT) : 1
				b[idx] += fac * b.w * exp(-0.5(cv - current_bin)^2 / b.bin_width^2)
			end	
		end

		return nothing
	end

	function update_bias_parametric!(b::BiasPotential, cv, itrj)
		batchsize = b.batchsize
		idx = mod1(itrj, batchsize)
		push!(b.cv_storage, cv)
		push!(b.bias_storage, b(cv))
		b.fullbias_storage[itrj] = b.values
		update_bias_regular!(b, cv)

		if is_symmetric(b)
			update_bias_regular!(b, -cv)
		end

		if idx == batchsize
			p = get_parameters(b)
			lb, ub = get_bounds(b)
			minimizer = get_minimizer(b)
			f = OptimizationFunction(b.testfun, Optimization.AutoForwardDiff())
			sol = solve(
				OptimizationProblem(f, p, b, lb = lb, ub = ub),
				minimizer,
				time_limit = 60,
			)
			b.current_parameters .= sol.u
			test = b.testfun(p, b)
			println_verbose(
				b.KS_fp,
				itrj, " ", sol.u[1], " ", sol.u[2], " ", sol.u[3], " ", test,
				"# itrj parameters test"
			)
			println("==============================================")
			flush(b.KS_fp)
		end

		return nothing
	end

	function parametricFES(parameters, cv)
		A, B, C = parameters
		return -A * cv^2 - B * cos(π * cv * C)^2
	end

	function parametricBias(parameters, cv)
		A, B, C = parameters
		return A * cv^2 + B * cos(π * cv * C)^2
	end

	function KStest(parameters, b::BiasPotential)
		#Sg = cv -> parametricFES(parameters, cv)
		batchsize = b.batchsize

		cv_data = get_cvstorage(b)
		len = length(cv_data)
		itvl = Int(len / batchsize)

		weight = (cv, bias, fullbias) -> weight_for_cdf(parameters, cv, bias, fullbias, b)
		bias_data = get_biasstorage(b)[1:itvl:end]
		fullbias_data = get_fullbiasstorage(b)[1:itvl:len]
		cv_data = cv_data[1:itvl:end]
		cdf = ecdf(cv_data, weights = weight.(cv_data, bias_data, fullbias_data))
		#cdf = ecdf(cv_data, weights = exp.(Sg.(cv_data) + bias_data))
		sorteddata = cdf.sorted_values

		cvmin = sorteddata[1]
		cvmax = sorteddata[end]

		l = abs(0.0 - (sorteddata[1] - cvmin) / (cvmax - cvmin))
		r = abs(cdf(sorteddata[1]) - (sorteddata[1] - cvmin) / (cvmax - cvmin))

		for i in 2:batchsize
			P_i = (sorteddata[i] - cvmin) / (cvmax - cvmin)
			lnext = abs(cdf(sorteddata[i-1]) - P_i)
			rnext = abs(cdf(sorteddata[i]) - P_i)
			l = l ≤ lnext ? lnext : l
			r = r ≤ rnext ? rnext : r
		end

		return max(l, r)
	end

	function GADtest(parameters, b::BiasPotential)		
		#Sg = cv -> parametricFES(parameters, cv)
		batchsize = b.batchsize

		cv_data = get_cvstorage(b)
		len = length(cv_data)
		itvl = Int(len / batchsize)

		weight = (cv, bias, fullbias) -> weight_for_cdf(parameters, cv, bias, fullbias, b)
		bias_data = get_biasstorage(b)[1:itvl:len]
		fullbias_data = get_fullbiasstorage(b)[1:itvl:len]
		cv_data = cv_data[1:itvl:end]

		#F = ecdf(cv_data, weights = exp.(Sg.(cv_data) + bias_data))
		F = ecdf(cv_data, weights = weight.(cv_data, bias_data, fullbias_data))
		sorteddata = F.sorted_values
		#cvmin, cvmax = get_CVlims(b)
		cvmin = floor(minimum(sorteddata))
		cvmax = ceil(maximum(sorteddata))

		S = 0.0
		
		for i in 1:batchsize - 1
			p_i = F(sorteddata[i])
        	u_i = wantedCDF(sorteddata[i], cvmin, cvmax)
        	u_ip1 = wantedCDF(sorteddata[i+1], cvmin, cvmax)
			S += p_i^2 * log(u_ip1 / u_i) - (p_i-1)^2 * log((1-u_ip1) / (1-u_i))
		end

		u_1 = wantedCDF(sorteddata[1], cvmin, cvmax)
		u_n = wantedCDF(sorteddata[end], cvmin, cvmax)
		return batchsize * (-1 - log(u_n) - log(1-u_1) + S)
	end

	function GADLTtest(parameters, b::BiasPotential)
		#Sg = cv -> parametricFES(parameters, cv)
		batchsize = b.batchsize

		cv_data = get_cvstorage(b)
		len = length(cv_data)
		itvl = Int(len / batchsize)

		weight = (cv, bias, fullbias) -> weight_for_cdf(parameters, cv, bias, fullbias, b)
		bias_data = get_biasstorage(b)[1:itvl:len]
		fullbias_data = get_fullbiasstorage(b)[1:itvl:len]
		cv_data = cv_data[1:itvl:end]

		#F = ecdf(cv_data, weights = exp.(Sg.(cv_data) + bias_data))
		F = ecdf(cv_data, weights = weight.(cv_data, bias_data, fullbias_data))
		sorteddata = F.sorted_values
		cvmin = floor(minimum(sorteddata))
		cvmax = ceil(maximum(sorteddata))

		S = 0.0

		for i in 1:batchsize - 1
			p_i = F(sorteddata[i])
			u_i = wantedCDF(sorteddata[i], cvmin, cvmax)
			u_ip1 = wantedCDF(sorteddata[i+1], cvmin, cvmax)
			S += p_i^2 * log(u_ip1/u_i) + 2p_i * (u_i - u_ip1)
		end

		p_n = F(sorteddata[end])
		u_n = wantedCDF(sorteddata[end], cvmin, cvmax)
		return batchsize * (0.5 - 2p_n * (1-u_n) - p_n^2 * log(u_n) + S)
	end

	function weight_for_cdf(parameters, cv, bias, fullbias, b::BiasPotential)
		cv_vals = b.cv_vals
		inorm = parametric_norm(parameters, fullbias, cv_vals)
		weight = inorm * exp( parametricFES(parameters, cv) + bias )
		return weight
	end
	
	function parametric_norm(parameters, old_bias, cv_vals)
		normA = zeros(nthreads() * 8)
		i = 0

		@batch for cv in cv_vals
			i += 1
			normA[threadid() * 8] += exp(-parametricFES(parameters, cv) - old_bias[i])
			#normA += exp( -parametricFES(parameters, cv) -
			#	parametricBias(old_bias, cv) )
		end

		return sum(normA)
	end

	function wantedCDF(cv::Float64, cvmin, cvmax)
		return (cv - cvmin) / (cvmax - cvmin)
	end
	
	function parametric_to_bias!(b::BiasPotential)
		parameters = b.current_parameters

		for (idx, cv) in enumerate(b.cv_vals)
			b[idx] = parametricBias(parameters, cv)
		end

		return nothing
	end

	function (b::BiasPotential)(cv::Float64)
		return return_potential(b, cv)
	end

	function return_potential(b::BiasPotential, cv)
		cvmin, cvmax = get_CVlims(b)

		if cvmin <= cv < cvmax
			grid_index = index(b, cv)
			return b[grid_index]
		else
			penalty = b.k * (0.1 + min((cv-cvmin)^2, (cv-cvmax)^2))
			return penalty
		end
	end

end


