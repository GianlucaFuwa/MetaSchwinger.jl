module Updates
    using LoopVectorization
    using Random

    import ..Gaugefields: Gaugefield, calc_Sg, recalc_CV!, substitute_U!, staple
    import ..Metadynamics: BiasPotential, update_bias!

    include("metro.jl")
    include("instanton.jl")
    include("tempering.jl")
end
