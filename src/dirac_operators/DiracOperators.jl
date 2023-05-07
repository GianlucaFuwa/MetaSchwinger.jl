module DiracOperatorModule

    using LinearAlgebra
    using Polyester
    using SparseArrays
    using StaticArrays

    import ..Gaugefields: Gaugefield
    import ..System_parameters: Params

    const eye2 = SMatrix{2,2,ComplexF64,4}([
        1 0
        0 1
    ])

    const γ1 = SMatrix{2,2,ComplexF64,4}([
        0 1
        1 0
    ])

    const γ2 = SMatrix{2,2,ComplexF64,4}([
        1 0
        0 -1
    ])

    function flat_index(NX, ix, it)
        idx = 1 + (ix - 1) * 2 + (it - 1) * NX * 2 
        return idx
    end
    
    function δ(i, j)
        return i == j
    end

    abstract type AbstractDiracOperator end

    include("naive_dirac.jl")
    include("wilson_dirac.jl")

    function dirac_operator(p::Params)
        if p.Dirac_operator === nothing 
            return nothing
        elseif p.Dirac_operator == "Naive"
            return NaiveDiracOperator(p.mass, p.BC)
        elseif p.Dirac_operator == "Wilson"
            return WilsonDiracOperator(p.mass, p.BC)
        else
            error("Dirac operator ", p.Dirac_operator, " not supported")
        end
    end
    
end
