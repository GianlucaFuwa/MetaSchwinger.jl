struct Kernel
    height::Float64
    center::Float64
    σ::Float64
end

# define evaluation of kernel both as functor and normal function to make it applicable
# to @turbo
(k::Kernel)(s, cutoff², penalty) = evaluate_kernel(s, k.height, k.center, k.σ, cutoff², penalty)

@inline function evaluate_kernel(s, height, center, σ, cutoff², penalty)
    diff = (center - s) / σ
    diff² = diff^2
    out = ifelse(diff²>=cutoff², 0.0, height * (exp(-0.5diff²) - penalty))
    return out
end

function Base.:*(c::Real, k::Kernel)
    return Kernel(c * k.height, k.center, k.σ)
end

function Base.:+(k::Kernel, other::Kernel) # overloaded addition for kernel merging
    h = k.height + other.height
    c = (k.height * k.center + other.height * other.center) / h
    s_my_part = k.height * (k.σ^2 + k.center^2)
    s_other_part = other.height * (other.σ^2 + other.center^2)
    s² = (s_my_part + s_other_part) / h - c^2
    return Kernel(h, c, sqrt(s²))
end

mutable struct OPES <: Bias
    explore::Bool
    is_first_step::Bool
    after_calculate::Bool

    is_static::Bool
    symmetric::Bool
    counter::Int64
    CVlims::NTuple{2,Float64}

    barrier::Float64
    biasfactor::Float64
    bias_prefactor::Float64
    stride::Int64

    σ₀::Float64
    adaptive_σ::Bool
    adaptive_σ_stride::Int64
    adaptive_counter::Int64
    av_cv::Union{Nothing, Float64}
    av_m2::Union{Nothing, Float64}
    σ_min::Float64
    fixed_σ::Bool

    ϵ::Float64
    sum_weights::Float64
    sum_weights²::Float64
    current_bias::Float64
    no_Z::Bool
    Z::Float64
    KDEnorm::Float64

    d_thresh::Float64
    cutoff²::Float64
    penalty::Float64
    nker::Int64
    kernels::Vector{Kernel}
    nδker::Int64
    δkernels::Vector{Kernel}
    kernelsfp::Union{Nothing, String}

    old_sum_weights::Float64
    old_Z::Float64
    old_KDEnorm::Float64

    statefp::Union{Nothing, IOStream}
    write_state_every::Int64

    function OPES(p::ParameterSet; instance=1)
        println(">> Setting OPES instance $(instance)...")
        is_first_step = true
        counter = 1
        Z = 1.0

        explore = false
        is_static = false
        symmetric = true

        # CVlims = p.CVlims===nothing ? (-Inf, Inf) : p.CVlims
        CVlims = (-6, 6)
        println("\t>> CVlims = $(CVlims)")
        @assert CVlims[1] < CVlims[2] "First entry in CVLIMS must be smaller than second"

        stride = p.stride
        println("\t>> STRIDE = $(stride)")
        barrier = p.barrier === nothing ? 0 : p.barrier
        @assert barrier >= 0 "BARRIER should be > 0"
        println("\t>> BARRIER = $(barrier)")
        biasfactor = p.biasfactor
        println("\t>> BIASFACTOR = $(biasfactor)")

        @assert biasfactor > 1 "BIASFACTOR must be > 1"
        bias_prefactor = 1 - 1 / biasfactor
        println("\t>> BIAS_PREFACTOR = $(bias_prefactor)")

        if explore
            @assert biasfactor != Inf "BIASFACTOR=Inf is not compatible with EXPLORE"
            bias_prefactor = biasfactor - 1
        end

        adaptive_σ = false
        adaptive_σ_stride = p.adaptive_σ_stride === nothing ? 0 : p.adaptive_σ_stride
        σ₀ = explore ? sqrt(barrier) : p.sigma0
        if σ₀ == 0.0
            adaptive_σ = true
            println("\t>> ADAPTIVE_SIGMA is used")
            @assert biasfactor != Inf "BIASFACTOR=Inf is not compatible with adaptive σ"
            adaptive_counter = 0
            adaptive_σ_stride == 0 && (adaptive_σ_stride = 100stride)
            av_cv = 0.0
            av_m2 = 0.0
            @assert adaptive_σ_stride >= stride "ADAPTIVE_SIGMA_STRIDE SHOULD BE >= STRIDE"
            println("\t>> ADAPTIVE_SIGMA_STRIDE = $(adaptive_σ_stride)")
        else
            @assert adaptive_σ_stride == 0 "if SIGMA0 is set then it cannot be adaptive"
            av_cv = nothing
            av_m2 = nothing
            adaptive_counter = 0
            if explore
                σ₀ = sqrt(biasfactor)
            end
            println("\t>> SIGMA0 = $(σ₀)")
        end

        σ_min = p.σ_min === nothing ? 1e-6 : p.σ_min
        println("\t>> MIN_SIGMA = $(σ_min)")

        ϵ = p.opes_epsilon === nothing ? exp(-barrier / bias_prefactor) : p.opes_epsilon
        @assert ϵ > 0 "EPSILON must be > 0. Your BARRIER might be too high"
        println("\t>> EPSILON = $(ϵ)")
        sum_weights = ϵ^bias_prefactor
        sum_weights² = sum_weights^2
        current_bias = 0.0

        cutoff = p.cutoff === nothing ? sqrt(2barrier / bias_prefactor) : p.cutoff
        @assert cutoff > 0 "CUTOFF must be > 0"
        println("\t>> CUTOFF = $(cutoff)")
        cutoff² = cutoff^2
        penalty = exp(-0.5cutoff²)

        nker = 0
        kernels = Vector{Kernel}(undef, 1000)
        nδker = 0
        δkernels = Vector{Kernel}(undef, 2)

        d_thresh = p.d_thresh === nothing ? 1.0 : p.d_thresh
        d_thresh != 0 && (@assert d_thresh > 0 && d_thresh^2 < cutoff² "")
        println("\t>> THRESHOLD = $(d_thresh)")

        no_Z = p.no_Z===nothing ? false : p.no_Z
        println("\t>> Z is $(ifelse(no_Z, "Static", "Dynamic"))")
        if no_Z
            sum_weights = 1.0
            sum_weights² = 1.0
        end

        fixed_σ = p.fixed_σ === nothing ? false : p.fixed_σ
        kernelsfp = p.kernelsfp === nothing ? nothing : p.kernelsfp

        statefp = p.statefp === nothing ? nothing : open(p.statefp, "w")
        write_state_every = p.write_state_every === nothing ? 0 : p.write_state_every
        if write_state_every != 0
            @assert statefp !== nothing "filename for storing states not specified"
        end
        if write_state_every > 0
            @assert write_state_every >= stride "write_state_every should be larger than stride"
        end

        old_sum_weights = sum_weights
        old_Z = Z
        KDEnorm = explore ? counter : sum_weights
        old_KDEnorm = KDEnorm

        if statefp !== nothing
            println(statefp, "sum_weights\tbias_prefactor\tZ\tepsilon\tnker\tsigma0")
            flush(statefp)
        end
        println()
        return new(
            explore, is_first_step, false,
            is_static, symmetric, counter, CVlims,
            barrier, biasfactor, bias_prefactor, stride,
            σ₀, adaptive_σ, adaptive_σ_stride, adaptive_counter, av_cv, av_m2, σ_min, fixed_σ,
            ϵ, sum_weights, sum_weights², current_bias, no_Z, Z, KDEnorm,
            d_thresh, cutoff², penalty, nker, kernels, nδker, δkernels, kernelsfp,
            old_sum_weights, old_Z, old_KDEnorm,
            statefp, write_state_every,
        )
    end
end

function (o::OPES)(cv) # OPES functor that calculates current_bias among other things
    if !in_bounds(cv, o.CVlims...)
        bounds_penalty = 1000
        bound = sign(cv)==1 ? o.CVlims[2] : o.CVlims[1]
        dist² = (cv - bound)^2
        calculate!(o, bound)
        return o.current_bias + bounds_penalty*dist²
    else
        calculate!(o, cv)
        return o.current_bias
    end
end

function calculate!(o::OPES, s)
    o.is_first_step && return nothing
    kernels = o.kernels

    prob = 0.0
    for i in 1:o.nker
        prob += kernels[i](s, o.cutoff², o.penalty)
    end
    prob /= o.sum_weights

    current_bias = o.bias_prefactor * log(prob / o.Z + o.ϵ)
    o.current_bias = current_bias
    o.after_calculate = true
    return nothing
end

function set_ϵ!(o::OPES, barrier)
    o.ϵ = exp(-barrier / o.bias_prefactor)
    return nothing
end

function set_cutoff!(o::OPES, barrier)
    o.cutoff² = 2barrier / o.bias_prefactor
    return nothing
end

function set_penalty!(o::OPES, val)
    o.penalty = val
    return nothing
end

function update_bias!(o::OPES, s; itrj=1)
    o.is_static && return nothing
    calculate!(o, s)
    update_opes!(o, s, itrj)
    if o.symmetric
        calculate!(o, -s)
        update_opes!(o, -s, itrj; symmetric=true)
    end
    return nothing
end

function update_opes!(o::OPES, s, itrj; symmetric=false)
    if o.is_first_step
        o.is_first_step = false
        return nothing
    end

    kernels = o.kernels
    δkernels = o.δkernels

    # update variance if adaptive σ
    if o.adaptive_σ && !symmetric
        o.adaptive_counter += 1
        τ = o.adaptive_σ_stride
        if o.adaptive_counter < o.adaptive_σ_stride
            τ = o.adaptive_counter
        end
        diff = s - o.av_cv
        o.av_cv += diff / τ
        o.av_m2 += diff * (s - o.av_cv)
        (o.adaptive_counter<o.adaptive_σ_stride && o.counter==1) && return nothing
    end

    if o.adaptive_σ && symmetric
        (o.adaptive_counter<o.adaptive_σ_stride && o.counter==1) && return nothing
    end

    (itrj%o.stride != 0 || !in_bounds(s, o.CVlims...)) && return nothing
    o.old_KDEnorm = o.KDEnorm
    old_nker = o.nker

    # get new kernel height
    height = exp(o.current_bias)

    # update sum_weights and neff
    o.counter += 1
    o.sum_weights += height
    o.sum_weights² += height^2
    neff = (1 + o.sum_weights)^2 / (1 + o.sum_weights²)
    if o.explore
        o.KDEnorm = o.counter
        height = 1.0
    else
        o.KDEnorm = o.sum_weights
    end

    # if needed rescale sigma and height
    σ = o.σ₀
    if o.adaptive_σ && !symmetric
        factor = o.explore ? 1 : o.biasfactor
        if o.counter == 2
            o.av_m2 *= o.biasfactor
            o.σ₀ = sqrt(o.av_m2 / o.adaptive_counter / factor)
            o.σ₀ = max(o.σ₀, o.σ_min)
        end
        σ = sqrt(o.av_m2 / o.adaptive_counter / factor)
        σ = max(σ, o.σ_min)
    end
    if !o.fixed_σ
        sz = o.explore ? o.counter : neff
        s_rescaling = (3sz / 4)^(-1/5)
        σ *= s_rescaling
        σ = max(σ, o.σ_min)
    end
    # height should be divided by sqrt(2π)*σ but this is cancelled out by Z
    # so we leave it out altogether but keep the s_rescaling
    height *= (o.σ₀ / σ)

    # add new kernel
    add_kernel!(o, height, s, σ, o.current_bias)

    # update Z
    if !o.no_Z
        # instead of redoing the full summation, we add only the changes, knowing that
        # uprob = old_uprob + delta_uprob
        # and we also need to consider that in the new sum there are some new centers
        # and some disappeared ones
        sum_uprob = 0.0
        δsum_uprob = 0.0
        for k in 1:o.nker
            for d in 1:o.nδker
                # take away contribution from kernels that are gone, and add new ones
                sgn = sign(δkernels[d].height)
                δsum_uprob += δkernels[d](kernels[k].center, o.cutoff², o.penalty) +
                    sgn * kernels[k](δkernels[d].center, o.cutoff², o.penalty)
            end
        end
        for d in 1:o.nδker
            for dd in 1:o.nδker
                # now subtract the δ_uprob added before, but not needed
                sgn = sign(δkernels[d].height)
                δsum_uprob -= sgn * δkernels[dd](δkernels[d].center, o.cutoff², o.penalty)
            end
        end

        sum_uprob = o.Z * o.old_KDEnorm * old_nker + δsum_uprob
        o.Z = sum_uprob / o.KDEnorm / o.nker
    end

    if o.write_state_every > 0 && itrj%o.write_state_every == 0 && !symmetric
        dump_state_to_file(o, itrj)
    end

    return nothing
end

function add_kernel!(o::OPES, height, s, σ, logweight=nothing)
    kernels = o.kernels
    δkernels = o.δkernels
    new_kernel = Kernel(height, s, σ)

    if o.d_thresh != 0
        i_min = get_mergeable_kernel(s, kernels, o.d_thresh, o.nker)

        if i_min < o.nker+1
            δkernels[1] = -1 * kernels[i_min]
            kernels[i_min] = kernels[i_min] + new_kernel
            δkernels[2] = kernels[i_min]
            o.nδker = 2
        else
            if o.nker+1 > length(kernels)
                push!(kernels, new_kernel)
                o.nker += 1
            else
                kernels[o.nker+1] = new_kernel
                o.nker += 1
            end
            δkernels[1] = new_kernel
            o.nδker = 1
        end
    end

    return nothing
end

function get_mergeable_kernel(s, kernels, d_thresh, nker)
    d_min = d_thresh
    i_min = nker+1

    for i in 1:nker
        d = abs(kernels[i].center-s) / kernels[i].σ
        ismin = d < d_min
        d_min = ifelse(ismin, d, d_min)
        i_min = ifelse(ismin, i, i_min)
    end

    return i_min
end

function dump_state_to_file(o::OPES, args...)
    (tmppath, tmpio) = mktemp()
    println(tmpio, "height\tcenter\tsigma")
    for kernel in o.kernels
        println(tmpio, "$(kernel.height)\t$(kernel.center)\t$(kernel.σ)")
    end
    close(tmpio)
    mv(tmppath, o.kernelsfp*".txt", force=true)

    println(o.statefp, "$(o.sum_weights)\t$(o.bias_prefactor)\t$(o.Z)\t$(o.ϵ)\t$(o.nker)\t$(o.σ₀)")
    flush(o.statefp)
    return nothing
end

function in_bounds(x::Real, lower_bound::Real, upper_bound::Real)
    lower_bound <= x <= upper_bound  && return true
    return false
end
