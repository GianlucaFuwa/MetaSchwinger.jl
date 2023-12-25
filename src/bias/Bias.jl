module BiasModule
    using Base.Threads
    using DelimitedFiles
    using ForwardDiff
    using Optimization
    using OptimizationOptimJL
    using StatsBase
    using StructArrays
    using Zygote

    import ..Parameters: ParameterSet
	import ..VerbosePrint: Verbose, Verbose_, println_verbose

    abstract type Bias end

    function update_bias!(::Bias, cv)
        error("No update! method for given Bias type")
    end

    function ∂V∂Q(::Bias, cv)
        error("No derivative method for given Bias type")
    end

    include("metadynamics.jl")
    include("opes.jl")

    ref(q) = -0.25q^2 + 19.5cos(1.045π*q)^2

    # function MSE(bias::Metadynamics)
    #     ref_bias = ref.(bias.cv_vals[2:end-1])
    #     ref_bias .-= minimum(ref_bias)
    #     current_bias = bias.values[2:end-1] .- minimum(bias.values[2:end-1])
    #     MSE = mean((ref_bias .- current_bias).^2)
    #     STD = std((ref_bias .- current_bias).^2)
    #     return MSE, STD
    # end

    function MSE(bias::Metadynamics)
        current_bias = bias.values[1:end-1] .- minimum(bias.values[1:end-1])
        MSE = mean(current_bias)
        STD = std(current_bias)
        return MSE, STD
    end

    # function MSE(bias::Metadynamics)
    #     δs = bias.bin_width
    #     ind5m = index(bias, -5.0)
    #     ind5p = index(bias, 5.0)
    #     ref_bias = ref.(-5:δs:5)
    #     ref_bias .-= minimum(ref_bias)
    #     current_bias = bias.values[ind5m:ind5p] .- minimum(bias.values)
    #     MSE = mean((ref_bias .- current_bias).^2)
    #     STD = std((ref_bias .- current_bias).^2)
    #     return MSE, STD
    # end

    function MSE(bias::OPES)
        cv_vals = range(bias.CVlims...; step=0.01)
        current_bias = bias.(cv_vals)
        current_bias .-= minimum(current_bias)
        MSE = mean(current_bias)
        STD = std(current_bias)
        return MSE, STD
    end

    # function MSE(bias::OPES)
    #     cv_vals = range(bias.CVlims...; step=0.01)
    #     ref_bias = ref.(cv_vals)
    #     ref_bias .-= minimum(ref_bias)
    #     current_bias = bias.(cv_vals)
    #     current_bias .-= minimum(current_bias)
    #     MSE = mean((ref_bias .- current_bias).^2)
    #     STD = std((ref_bias .- current_bias).^2)
    #     return MSE, STD
    # end
end
