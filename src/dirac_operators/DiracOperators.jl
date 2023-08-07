module DiracOperators
    using LinearAlgebra
    using Polyester
    using SparseArrays
    using StaticArrays

    import ..Gaugefields: Gaugefield
    import ..Parameters: ParameterSet

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

    function dirac_operator(p::ParameterSet)
        if p.operator === nothing 
            return nothing
        elseif p.operator == "Naive"
            return NaiveDiracOperator(p.N, p.mass, p.BC)
        elseif p.operator == "Wilson"
            return WilsonDiracOperator(p.N, p.mass, p.BC)
        else
            error("Dirac operator ", p.operator, " not supported")
        end
    end

    function Base.setindex!(D::AbstractDiracOperator, v, i1, i2)
		@inbounds D.Dop[i1,i2] = v
	end

    function Base.setindex!(D::AbstractDiracOperator, v, i)
		@inbounds D.Dop[i] = v
	end

	@inline function Base.getindex(D::AbstractDiracOperator, i1, i2)
		@inbounds return D.Dop[i1,i2]
	end

    @inline function Base.getindex(D::AbstractDiracOperator, i)
		@inbounds return D.Dop[i]
	end

    function eachindex(D::AbstractDiracOperator)
        return eachindex(D.Dop)
    end

	function Base.size(D::AbstractDiracOperator)
		return size(D.Dop)
	end

    function clear!(D::AbstractDiracOperator)
        @batch for ii in eachindex(D)
            D[ii] = 0.0+0.0im
        end
        return nothing
    end
    
end
