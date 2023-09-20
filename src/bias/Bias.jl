module BiasModule
    using Base.Threads
    using DelimitedFiles
    using LoopVectorization
    using ForwardDiff
    using Optimization
    using OptimizationOptimJL
    using Polyester
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

end
