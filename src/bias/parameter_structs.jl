Base.@kwdef mutable struct MetaDParameters{O}
    is_static::Bool = false
    symmetric::Bool = true
    CVlims::NTuple{2, Float64} = (-7, 7)
    bin_width::Float64 = 1e-2
    w::Float64 = 1e-3
    k::Float64 = 1000
    parametric::Bool = false
    current_parameters::Union{Nothing, Vector{Float64}} = nothing
    lower_bounds::Union{Nothing, Vector{Float64}} = nothing
    upper_bounds::Union{Nothing, Vector{Float64}} = nothing
    batchsize::Union{Nothing, Int64} = nothing
    cv_storage::Union{Nothing, Vector{Float64}} = nothing
    bias_storage::Union{Nothing, Vector{Float64}} = nothing
    fullbias_storage::Union{Nothing, Vector{Vector{Float64}}} = nothing
    testfun::Union{Nothing, Function} = nothing
    minimizer::O = nothing
    ΔT::Float64 = Inf
    values::Vector{Float64} = zeros(1 + round(Int, (CVlims[2] - CVlims[1]) / bin_width))
    cv_vals::Vector{Float64} = range(CVlims[1], CVlims[2], step = bin_width)
    exceeded_count::Int64 = 0
    fp::Union{Nothing, IOStream} = nothing
    KS_fp::Union{Nothing, Verbose} = nothing
end

Base.@kwdef mutable struct OPESParameters
    explore::Bool = false
    is_first_step::Bool = true
    after_calculate::Bool = false

    is_static::Bool = false
    symmetric::Bool = false
    counter::Int64 = 1

    barrier::Float64 = 40
    biasfactor::Float64 = barrier
    bias_prefactor::Float64 = 1-1/biasfactor
    stride::Int64 = 1

    σ₀::Float64 = 0.01
    adaptive_σ::Bool = false
    adaptive_σ_stride::Int64 = 0
    adaptive_counter::Int64 = 0
    av_cv::Float64 = 0.0
    av_m2::Float64 = 0.0
    σ_min::Float64 = 1e-6
    fixed_σ::Bool = false

    ϵ::Float64 = exp(-barrier/bias_prefactor)
    sum_weights::Float64 = ϵ^bias_prefactor
    sum_weights²::Float64 = sum_weights^2
    current_bias::Float64 = 0.0
    no_Z::Bool = false
    Z::Float64 = 1
    KDEnorm::Float64 = sum_weights

    d_thresh::Float64 = 1.0
    cutoff²::Float64 = 2barrier/bias_prefactor
    penalty::Float64 = exp(-0.5cutoff²)
    kernels::StructArray{Kernel} = StructArray(Vector{Kernel}(undef, 0))
    kernelsfp::Union{Nothing, IOStream} = nothing

    old_sum_weights::Float64 = sum_weights
    old_Z::Float64 = Z
    old_KDEnorm::Float64 = KDEnorm
    δkernels::StructArray{Kernel} = StructArray(Vector{Kernel}(undef, 0))

    probfp::Union{Nothing, IOStream} = nothing
    write_prob_stride::Int64 = 1
    store_old_probs::Bool = false
end

Base.@kwdef mutable struct OPESParameters1
    explore::Bool = true
    is_first_step::Bool = true
    after_calculate::Bool = false

    is_static::Bool = false
    symmetric::Bool = false
    counter::Int64 = 1

    barrier::Float64 = 40.0
    bias_prefactor::Float64 = barrier-1
    stride::Int64 = 1

    σ₀::Float64 = 0.1*sqrt(barrier)
    adaptive_σ::Bool = false
    adaptive_σ_stride::Int64 = 0
    adaptive_counter::Int64 = 0
    av_cv::Float64 = 0.0
    av_m2::Float64 = 0.0
    σ_min::Float64 = 1e-6
    fixed_σ::Bool = false

    ϵ::Float64 = exp(-barrier/bias_prefactor)
    sum_weights::Float64 = ϵ^bias_prefactor
    sum_weights²::Float64 = sum_weights^2
    current_bias::Float64 = 0.0
    no_Z::Bool = false
    Z::Float64 = 1
    KDEnorm::Float64 = counter

    d_thresh::Float64 = 1.0
    cutoff²::Float64 = 2*barrier
    penalty::Float64 = exp(-0.5*cutoff²)
    kernels::StructArray{Kernel} = StructArray(Vector{Kernel}(undef, 0))
    kernelsfp::Union{Nothing, String} = nothing

    old_sum_weights::Float64 = sum_weights
    old_Z::Float64 = Z
    old_KDEnorm::Float64 = KDEnorm
    δkernels::StructArray{Kernel} = StructArray(Vector{Kernel}(undef, 0))

    probfp::Union{Nothing, IOStream} = nothing
    write_prob_stride::Int64 = 1
    store_old_probs::Bool = false
end

function set_params_value!(value_Params, values, pnames)
    d = struct2dict(values)

    for (i, pname_i) in enumerate(pnames)
        if haskey(d,String(pname_i))
            if d[String(pname_i)] == "nothing"
                value_Params[i] = nothing
            else
                value_Params[i] = d[String(pname_i)]
            end
        end
    end

    return nothing
end

function struct2dict(x::T) where {T}
    dict = Dict{String,Any}(string(fn) => getfield(x, fn) for fn in fieldnames(T))
    return dict
end
