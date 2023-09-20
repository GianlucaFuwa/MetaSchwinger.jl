struct Metadynamics{O} <: Bias
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
    biasfile::Union{Nothing, String}
    write_state_every::Union{Nothing, Int64}
    KS_fp::Union{Nothing, Verbose}

    function Metadynamics(p::ParameterSet; instance=1)
        println(">> Setting MetaD instance $(instance)...")
        instance === 0 && println("\t>> instance 0 has dummy bias")
        is_static = instance==0 ? true : p.is_static[instance]
        symmetric = p.symmetric

        CVlims = p.CVlims
        instance !== 0 && println("\t>> CVLIMS = $(CVlims)")
        bin_width = p.bin_width
        instance !== 0 && println("\t>> BIN_WIDTH = $(bin_width)")
        w = p.w
        instance !== 0 && println("\t>> WEIGHT = $(w)")
        k = p.k
        instance !== 0 && println("\t>> PENALTY_WEIGHT = $(k)")

        parametric = p.parametric
        parametric && println("\t>> Setting Parametric MetaD")

        if parametric == true && is_static == false
            current_parameters = p.potential_parameters
            lower_bounds = p.lower_bounds
            upper_bounds = p.upper_bounds
            batchsize = p.batchsize
            cv_storage = Float64[]
            bias_storage = Float64[]
            fullbias_storage = Vector{Vector{Float64}}(undef, 0)
            testfun = getfield(BiasModule, Symbol(p.testfun))
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
        ΔT = p.ΔT===nothing ? Inf : p.ΔT
        instance !== 0 && println("\t>> BIASFACTOR = $(ΔT)")

        if instance == 0
            values = potential_from_file(p, nothing)
            instance !== 0 && println("\t>> BIAS INITIATED AS ZEROS")
        else
            values = potential_from_file(p, p.usebiases[instance])
        end

        cv_vals = range(CVlims[1], CVlims[2], step=bin_width)

        exceeded_count = 0
        biasfile = instance==0 ? nothing : p.biasfiles[instance]
        write_state_every = instance==0 ? nothing : p.write_state_every

        if p.parametric == true
            KS_fp = Verbose_(open(
                pwd() * p.measure_dir * "/" * p.testfun * "_" * p.minimizer *".txt",
                "w",
            ))
        else
            KS_fp = nothing
        end
        println()
        return new{typeof(minimizer)}(
            is_static, symmetric,
            CVlims, bin_width, w, k,
            parametric, current_parameters, lower_bounds, upper_bounds, batchsize,
            cv_storage, bias_storage, fullbias_storage, testfun, minimizer,
            well_tempered, ΔT,
            values, cv_vals,
            exceeded_count, biasfile, write_state_every, KS_fp,
        )
    end
end

include("parametric_metad.jl")

function potential_from_file(p::ParameterSet, usebias::Union{Nothing,String})
    if usebias === nothing
        len = 1 + round(
            Int,
            (p.CVlims[2] - p.CVlims[1]) / p.bin_width,
            RoundNearestTiesAway
        )
        println("\t>> BIAS INITIATED AS ZEROS")
        return zeros(len)
    else
        values = readdlm(usebias, Float64, comments = true)
        len = 1 + round(
            Int,
            (p.CVlims[2] - p.CVlims[1]) / p.bin_width,
            RoundNearestTiesAway
        )
        @assert length(values[:,2]) == len  "Potential length doesn't match parameters"
        println("\t>> BIAS INITIATED FROM: \"$(usebias)\"")
        return values[:,2]
    end
end

function Base.setindex!(b::Metadynamics, v, i)
    b.values[i] = v
    return nothing
end

@inline function Base.getindex(b::Metadynamics, i)
    return b.values[i]
end

@inline function index(b::Metadynamics, cv)
    grid_index = (cv - b.CVlims[1]) / b.bin_width + 0.5
    return round(Int, grid_index, RoundNearestTiesAway)
end

function Base.sum(b::Metadynamics)
    return sum(b.values)
end

function clear!(b::Metadynamics)
    for i in eachindex(b.values)
        b.values[i] = 0
    end

    return nothing
end

function update_bias!(b::Metadynamics, cv; itrj=1)
    b.is_static && return nothing
    if b.parametric == false
        update_bias_regular!(b, cv)

        if b.symmetric
            update_bias_regular!(b, -cv)
        end

    elseif b.parametric == true
        update_bias_parametric!(b, cv, itrj)
    end

    if b.write_state_every>0 && itrj%b.write_state_every == 0
        dump_state_to_file(b, itrj)
    end
    return nothing
end

function update_bias_regular!(b::Metadynamics, cv)
    grid_index = index(b, cv)

    if 1 <= grid_index < length(b.values)
        for (idx, current_bin) in enumerate(b.cv_vals)
            fac = b.well_tempered ? exp(-b[idx] / b.ΔT) : 1
            b[idx] += fac * b.w * exp(-0.5(cv - current_bin)^2 / b.bin_width^2)
        end
    end

    return nothing
end

function (b::Metadynamics)(cv::Float64)
    return return_bias(b, cv)
end

function return_bias(b::Metadynamics, cv)
    cvmin, cvmax = b.CVlims

    if cvmin <= cv < cvmax
        grid_index = index(b, cv)
        return b[grid_index]
    else
        penalty = b.k * (0.1 + min((cv - cvmin)^2, (cv - cvmax)^2))
        return penalty
    end
end

function ∂V∂Q(b::Metadynamics, cv)
    bw = b.bin_width
    num = -b(cv+2bw) + 8b(cv+bw) - 8b(cv-bw) + b(cv-2bw)
    denom = 12bw
    return num / denom
end

function write_to_file(b::Metadynamics, fp)
    for cv in b.cv_vals
        value = b(cv)
        println(fp, cv, "\t", value)
    end
    return nothing
end

function dump_state_to_file(b::Metadynamics, itrj)
    (tmppath, tmpio) = mktemp()
    println(tmpio, "cv\tV(cv)")

    for cv in b.cv_vals
        value = b(cv)
        println(tmpio, "$(cv)\t$(value)")
    end

    close(tmpio)
    mv(tmppath, b.biasfile*"_$itrj.txt", force = true)
    return nothing
end
