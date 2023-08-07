module OPESModule
	using DelimitedFiles
    using LoopVectorization
	using Polyester
	using StructArrays

	export OPES, Kernel, update!, calculate!

	struct ParameterSet end
	struct Verbose end
	# import ..Parameters: ParameterSet
	# import ..VerbosePrint: Verbose, Verbose_, println_verbose

	struct Kernel
		height::Float64
		center::Float64
		σ::Float64
		cutoff²::Float64
		penalty::Float64

		function Kernel(h, center, σ, cutoff², penalty)
			return new(h, center, σ, cutoff², penalty)
		end
	end

	@inline function (k::Kernel)(s) # Kernel functor that calculates Gaussian
		diff = (k.center - s) / k.σ
		diff² = diff^2
		diff² >= k.cutoff² && return 0.0
		return k.height * (exp(-0.5diff²) - k.penalty)
	end

	function Base.:*(c::Real, k::Kernel)
		return Kernel(c * k.height, k.center, k.σ, k.cutoff², k.penalty)
	end

	function Base.:+(k::Kernel, other::Kernel) # overloaded addition for kernel merging
		h = k.height + other.height
		c = (k.height * k.center + other.height * other.center) / h
		s_my_part = k.height * (k.σ^2 + k.center^2)
		s_other_part = other.height * (other.σ^2 + other.center^2)
		s² = (s_my_part + s_other_part) / h - c^2
		return Kernel(h, c, sqrt(s²), k.cutoff², k.penalty)
	end

    mutable struct OPES
        is_static::Bool
		symmetric::Bool
		counter::Int64

		bias_prefactor::Float64
		stride::Int64
		σ₀::Float64
		adaptive_σ_stride::Int64
		adaptive_counter::Int64
		fixed_σ::Bool
		ϵ::Float64
		sum_weights::Float64
		sum_weights²::Float64
		current_bias::Float64
		no_Z::Bool
		Z::Float64

		d_thresh::Float64
		cutoff²::Float64
		penalty::Float64
		kernels::StructArray{Kernel}
		kernelsfp::Union{Nothing, IOStream}

		work::Float64
		old_sum_weights::Float64
		old_Z::Float64
		δkernels::StructArray{Kernel} # keep track of the change in kernels after every step

		probfp::Union{Nothing, IOStream}
		write_prob_stride::Int64
		store_old_probs::Bool

		function OPES(; p::Union{Nothing, ParameterSet} = nothing)
			if p === nothing
				pnames = fieldnames(OPES)
				numparams = length(pnames)
				value_Params = Vector{Any}(undef, numparams)
				opes = OPESParameters()
				set_params_value!(value_Params, opes, pnames)
				return new(value_Params...)
			end

			counter = 1
			Z = 1.0
			work = 0.0
			is_static = p.is_static
			symmetric = p.symmetric

			stride = p.stride
			barrier = p.barrier; @assert barrier >= 0 "BARRIER should be > 0"
			biasfactor = p.biasfactor; @assert biasfactor > 1 "BIASFACTOR must be > 1"
			bias_prefactor = 1 - 1 / biasfactor
			adaptive_σ_stride = p.adaptive_σ_stride === nothing ? 0 : p.adaptive_σ_stride
			adaptive_σ_stride == 0 && (adaptive_σ_stride = 10stride)
			σ₀ = p.sigma0
			adaptive_counter = 0


			ϵ = p.opes_epsilon === nothing ? exp(-barrier / bias_prefactor) : p.opes_epsilon
			@assert ϵ > 0 "EPSILON must be > 0. Your BARRIER might be too high"
			sum_weights = ϵ^bias_prefactor
			sum_weights² = sum_weights^2

			cutoff = p.cutoff === nothing ? sqrt(2barrier / bias_prefactor) : p.cutoff
			@assert cutoff > 0 "CUTOFF must be > 0"
			cutoff² = cutoff^2
			penalty = exp(-0.5cutoff²)

			d_thresh = p.d_thresh === nothing ? 1.0 : p.d_thresh
			d_thresh != 0 && (@assert d_thresh > 0 && d_thresh^2 < cutoff² "")

			no_Z = p.no_Z === nothing ? false : p.no_Z
			if no_Z
				sum_weights = 1.0
				sum_weights² = 1.0
			end

			fixed_σ = p.fixed_σ === nothing ? false : p.fixed_σ
			kernelsfp = p.kernelsfp === nothing ? nothing : open(p.kernelsfp, "w")

			probfp = p.probfp === nothing ? nothing : open(p.probfp, "w")
			store_old_probs = p.store_old_probs === nothing ? true : p.store_old_probs
			write_prob_stride = p.write_prob_stride === nothing ? 0 : p.write_prob_stride

			old_sum_weights = sum_weights
			old_Z = Z
			kernels = StructArray(Vector{Kernel}(undef, 0))
			δkernels = StructArray(Vector{Kernel}(undef, 0))

			return new(
				is_static, symmetric, counter,
				bias_prefactor, stride, σ₀, adaptive_σ_stride, adaptive_counter, fixed_σ, ϵ,
				sum_weights, sum_weights², current_bias, Z,
				d_thresh, cutoff², penalty, kernels, kernelsfp,
				work, old_sum_weights, old_Z, δkernels,
				probfp, write_prob_stride, store_old_probs,
			)
		end
    end

	function (o::OPES)(cv) # OPES functor that calculates current_bias among other things
		calculate!(o, cv)
		return o.current_bias
	end

	function calculate!(o::OPES, cv)
		prob = 0.0
		@batch for kernel in o.kernels
			prob += kernel(cv)
		end

		prob /= o.sum_weights
		current_bias = o.bias_prefactor * log(prob / o.Z + o.ϵ)
		o.current_bias = current_bias

		# calculate work
		tot_δ = 0
		@batch for δkernel in o.δkernels
			tot_δ += δkernel(cv)
		end

		old_prob = (prob * o.sum_weights - tot_δ) / o.old_sum_weights
		o.work += current_bias - o.bias_prefactor * log(old_prob / o.old_Z + o.ϵ)
		return nothing
	end

	function update!(o::OPES, s; itrj = 1)
		kernels = o.kernels
		# update variance if adaptive σ

		# work done by the bias in one iteration, sues as zero reference a point at Inf
		# so that the work is always positive
		min_shift = o.bias_prefactor * log(o.old_Z/o.Z * o.old_sum_weights/o.sum_weights)
		o.work = 0
		o.δkernels = StructArray(Vector{Kernel}(undef, 0))
		o.old_sum_weights = o.sum_weights
		o.old_Z = o.Z
		old_nker = length(kernels)

		δkernels = o.δkernels
		# get new kernel height
		height = exp(o.current_bias)

		# update sum_weights and neff
		sum_heights = height
		sum_heights² = height^2
		o.counter += 1
		o.sum_weights += sum_heights
		o.sum_weights² += sum_heights²
		neff = (1 + o.sum_weights)^2 / (1 + sum_heights²)

		# if needed rescale sigma and height
		σ = o.σ₀

		if !o.fixed_σ
			s_rescaling = (3neff / 4)^(-1/5)
			σ = s_rescaling
			# height should be divided by sqrt(2π)*σ but this is canelled out by Z
			# so we leave it out altogether but keep the s_rescaling
			height /= s_rescaling
		end
		# add new kernel
		add_kernel!(o, height, s, σ, true)
		# update Z
		if !o.no_Z
			sum_uprob = 0
			δsum_uprob = 0
			@batch for k in eachindex(kernels)
				for d in eachindex(δkernels)
					# take away contribution from kernels that are gone, and add new ones
					sgn = δkernels[d].height < 0 ? -1 : 1
					δsum_uprob += δkernels[d](kernels[k].center)
					δsum_uprob += sgn * kernels[k](δkernels[d].center)
				end
			end
			sum_uprob = o.Z * o.old_sum_weights * old_nker + δsum_uprob
			o.Z = sum_uprob / sum_weights / length(kernels)
		end

		return nothing
	end

	function add_kernel!(o::OPES, height, s, σ, write_to_file)
		kernels = o.kernels
		δkernels = o.δkernels
		new_kernel = Kernel(height, s, σ, o.cutoff², o.penalty)

		if o.d_thresh != 0
			i_min, d_min = get_mergeable_kernel(s, kernels)

			if d_min < o.d_thresh
				push!(δkernels, -1 * kernels[i_min])
				kernels[i_min] = kernels[i_min] + new_kernel
				push!(δkernels, kernels[i_min])
			else
				push!(kernels, new_kernel)
				push!(δkernels, new_kernel)
			end
		end

		if write_to_file && o.kernelsfp !== nothing
			println(o.kernelsfp, "$(height)\t$(center)\t$(σ)\t$(o.current_bias)")
		end

		return nothing
	end

	function get_mergeable_kernel(s, kernels)
		d_min = Inf
		i_min = 0
		centers = view(kernels.center, :)
		σs = view(kernels.σ, :)

		@turbo for i in eachindex(kernels)
			d = abs(centers[i] - s) / σs[i]
			ismin = d < d_min
			d_min = ifelse(ismin, d, d_min)
			i_min = ifelse(ismin, i, i_min)
		end

		return i_min, d_min
	end

	include("parameter_structs.jl")

end
