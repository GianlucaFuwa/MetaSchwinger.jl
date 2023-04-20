module Metadynamics
	using DelimitedFiles
	using Optimization, ForwardDiff, Zygote, OptimizationOptimJL, StatsBase
	import ..System_parameters: Params
	import ..Verbose_print: Verbose, Verbose_, println_verbose
	import ..Gaugefields: Gaugefield
	
	struct Bias_potential{O}
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
		
		function Bias_potential(p::Params; instance::Int64 = 1)
			is_static = instance==0 ? true : p.is_static[instance]
			symmetric = p.symmetric

			CVlims = p.CVlims
			bin_width = p.bin_width
			w = p.w
			k = p.k

			parametric = p.parametric
			if parametric==true && is_static==false
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

			values = instance==0 ? potential_from_file(p,nothing) : potential_from_file(p, p.usebiases[instance])
			cv_vals = range(CVlims[1], CVlims[2], step=bin_width)

			exceeded_count = 0
			fp = instance==0 ? nothing : open(p.biasfiles[instance], "w")
			KS_fp = p.parametric==true ? Verbose_(open(pwd()*p.savebias_dir*"/"*p.testfun*"_tests_"*p.minimizer*".txt", "w")) : nothing

			return new{typeof(minimizer)}(
			is_static, symmetric,
			CVlims, bin_width, w, k,
			parametric, current_parameters, lower_bounds, upper_bounds, batchsize, cv_storage, bias_storage, fullbias_storage, testfun, minimizer,
			well_tempered, ΔT,
			values, cv_vals,
			exceeded_count, fp, KS_fp,
			)
		end	
	end

	function potential_from_file(p::Params,usebias::Union{Nothing,String})
		if usebias === nothing
			return zeros(round(Int,(p.CVlims[2]-p.CVlims[1])/p.bin_width,RoundNearestTiesAway)+1)
		else
			values = readdlm(usebias, Float64, comments=true)
			len = round(Int,(p.CVlims[2]-p.CVlims[1])/p.bin_width,RoundNearestTiesAway)+1
			@assert length(values[:,2]) == len  "Potential length should be $len but is $(length(values[:,2]))"
			return values[:,2]
		end
	end

	function write_to_file!(b::Bias_potential)
		for (idx, cv) in enumerate(b.cv_vals)
			value = ReturnPotential(b,cv)
			println(b.fp, "$cv $value # Metapotential")
		end
		return nothing
	end

	function Base.flush(b::Bias_potential)
        if b.fp !== nothing
            flush(b.fp)
        end
    end

	function Base.seekstart(b::Bias_potential)
		if b.fp !== nothing
            seekstart(b.fp)
        end
    end

	function Base.setindex!(b::Bias_potential,v,i::Int)
		b.values[min(length(b.values), max(i, 1))] = v
		return nothing
	end

	@inline function Base.getindex(b::Bias_potential,i::Int)
		return b.values[min(length(b.values), max(i, 1))]
	end

	@inline function index(b::Bias_potential,cv::Float64)
		grid_index = (cv - b.CVlims[1]) / b.bin_width + 0.5
		return round(Int,grid_index,RoundNearestTiesAway)
	end

	function Base.sum(b::Bias_potential)
		return sum(b.values)
	end

	is_static(b::Bias_potential) = b.is_static
	is_symmetric(b::Bias_potential) = b.symmetric
	is_well_tempered(b::Bias_potential) = b.well_tempered
	is_parametric(b::Bias_potential) = b.parametric

	get_CVlims(b::Bias_potential) = b.CVlims
	get_parameters(b::Bias_potential) = b.current_parameters
	get_bounds(b::Bias_potential) = b.lower_bounds, b.upper_bounds
	get_minimizer(b::Bias_potential) = b.minimizer
	get_batchsize(b::Bias_potential) = b.batchsize
	get_cvstorage(b::Bias_potential) = b.cv_storage
	get_biasstorage(b::Bias_potential) = b.bias_storage
	get_fullbiasstorage(b::Bias_potential) = b.fullbias_storage

	function update_bias!(b::Bias_potential, cv::Float64; itrj=nothing)
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

	function update_bias_regular!(b::Bias_potential, cv::Float64)
		grid_index = index(b, cv)
		if 1 <= grid_index < length(b.values)
			for (idx, current_bin) in enumerate(b.cv_vals)
				well_tempered_fac = is_well_tempered(b) ? exp(-b[idx] / b.ΔT) : 1
				b[idx] += well_tempered_fac*b.w*exp(-0.5(cv-current_bin)^2 / b.bin_width^2)
			end	
		end
		return nothing
	end

	function update_bias_parametric!(b::Bias_potential, cv::Float64 ,itrj)
		batchsize = get_batchsize(b)
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
			sol = solve(OptimizationProblem(f, p, b, lb=lb, ub=ub), minimizer, time_limit=60)
			b.current_parameters .= sol.u
			test = b.testfun(p, b)
			println_verbose(b.KS_fp, "$itrj $(sol.u[1]) $(sol.u[2]) $(sol.u[3]) $test # itrj parameters test")
			println("========================================")
			flush(b.KS_fp)
		end
		return nothing
	end

	function parametricFES(parameters, cv::Float64)
		A, B, C = parameters
		return -A*cv^2 - B*cos(pi*cv*C)^2
	end

	function parametricBias(parameters, cv::Float64)
		A, B, C = parameters
		return A*cv^2 + B*cos(pi*cv*C)^2
	end

	function KStest(parameters, b::Bias_potential)
		Sg = cv -> parametricFES(parameters, cv)
		batchsize = get_batchsize(b)

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

		l = abs( 0.0 - (sorteddata[1] - cvmin) / (cvmax - cvmin) )
		r = abs( cdf(sorteddata[1]) - (sorteddata[1] - cvmin) / (cvmax - cvmin) )
		for i = 2:batchsize
			P_i = (sorteddata[i] - cvmin) / (cvmax - cvmin)
			lnext = abs(cdf(sorteddata[i-1]) - P_i)
			rnext = abs(cdf(sorteddata[i]) - P_i)
			l = l ≤ lnext ? lnext : l
			r = r ≤ rnext ? rnext : r
		end
		return max(l, r)
	end

	function GADtest(parameters, b::Bias_potential)		
		Sg = cv -> parametricFES(parameters, cv)
		batchsize = get_batchsize(b)

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
		for i = 1:batchsize-1
			p_i = F(sorteddata[i])
        	u_i = wantedCDF(sorteddata[i], cvmin, cvmax)
        	u_ip1 = wantedCDF(sorteddata[i+1], cvmin, cvmax)
			S += p_i^2 * log(u_ip1 / u_i) -
            (p_i-1)^2 * log((1-u_ip1) / (1-u_i))
		end
		u_1 = wantedCDF(sorteddata[1], cvmin, cvmax)
		u_n = wantedCDF(sorteddata[end], cvmin, cvmax)
		return batchsize * ( -1 - log(u_n) - log(1-u_1) + S )
	end

	function GADLTtest(parameters, b::Bias_potential)
		Sg = cv -> parametricFES(parameters, cv)
		batchsize = get_batchsize(b)

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
		for i = 1:batchsize-1
			p_i = F(sorteddata[i])
			u_i = wantedCDF(sorteddata[i], cvmin, cvmax)
			u_ip1 = wantedCDF(sorteddata[i+1], cvmin, cvmax)
			S += p_i^2 * log(u_ip1/u_i) + 2p_i * (u_i - u_ip1)
		end
		p_n = F(sorteddata[end])
		u_n = wantedCDF(sorteddata[end], cvmin, cvmax)
		return batchsize * ( 0.5 - 2p_n*(1-u_n) - p_n^2*log(u_n) + S )
	end

	function weight_for_cdf(parameters, cv, bias, fullbias, b::Bias_potential)
		cv_vals = b.cv_vals

		inorm = parametric_norm(parameters, fullbias, cv_vals)
		weight = inorm * exp( parametricFES(parameters, cv) + bias )
		return weight
	end
	
	function parametric_norm(parameters, old_bias, cv_vals)
		normA = 0.0
		i = 0
		for cv in cv_vals
			i += 1
			normA += exp( -parametricFES(parameters, cv) -old_bias[i] )
			#normA += exp( -parametricFES(parameters, cv) -
			#	parametricBias(old_bias, cv) )
		end
		return normA
	end

	function wantedCDF(cv::Float64, cvmin, cvmax)
		return (cv - cvmin) / (cvmax - cvmin)
	end
	
	function parametric_to_bias!(b::Bias_potential)
		parameters = get_parameters(b)
		for (idx, cv) in enumerate(b.cv_vals)
			b[idx] = parametricBias(parameters, cv)
		end
		return nothing
	end

	function (b::Bias_potential)(cv::Float64)
		return ReturnPotential(b, cv)
	end

	function ReturnPotential(b::Bias_potential, cv::Float64)
		cvmin, cvmax = get_CVlims(b)
		if cvmin <= cv < cvmax
			grid_index = index(b, cv)
			return b[grid_index]
		else
			penalty = b.k * ( 0.1 + min( (cv-cvmin)^2, (cv-cvmax)^2 ) )
			return penalty
		end
	end

	function DeltaV(b::Bias_potential, cvold::Float64, cvnew::Float64)
		dV = b(cvnew) - b(cvold)
		return dV
	end

end


