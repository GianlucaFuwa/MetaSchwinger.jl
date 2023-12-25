module Gaugefields
	using LoopVectorization
	import ..Parameters: ParameterSet

    abstract type Abstractfield end

	mutable struct Gaugefield <: Abstractfield
		NX::Int64
		NT::Int64
		NV::Int64
		β::Float64
		U::Array{Float64,3}
		Sg::Float64
		CV::Float64

		function Gaugefield(p::ParameterSet; instance=1)
            println(">> Setting Gaugefield instance $(instance)...")
			NX = p.N[1]
			NT = p.N[2]
            println("\t>> L = $((NX, NT))")
			NV = NX * NT
			β = p.β
            println("\t>> beta = $(β)")
			U = zeros(Float64, 2, NX, NT)
            println()
			return new(NX, NT, NV, β, U, 0.0, 0.0)
		end

		function Gaugefield(NX, NT, β)
			NV = NX * NT
			U = zeros(Float64, 2, NX, NT)
			return new(NX, NT, NV, β, U, 0.0, 0.0)
		end
	end

	function Base.setindex!(u::Abstractfield, v, μ, ix, it)
		@inbounds u.U[μ,ix,it] = v
	end

	@inline function Base.getindex(u::Abstractfield, μ, ix, it)
		@inbounds return u.U[μ,ix,it]
	end

	function Base.size(u::Abstractfield)
		return u.NX, u.NT
	end

    function Base.similar(u::T) where {T<:Abstractfield}
        T==Liefield && return T(size(u)...)
        T==Gaugefield && return T(size(u)..., u.β)
    end

	function Base.cis(u::Gaugefield)
		NX, NT = size(u)
		expiU = zeros(ComplexF64, 2, NX, NT)

        for it in 1:NT
            for ix in 1:NX
                for μ in 1:2
					expiU[μ,ix,it] = cis(u[μ,ix,it])
				end
			end
		end

		return expiU
	end

	function Base.circshift(u::Gaugefield, shifts::NTuple{2,Int}, μ)
		return circshift(u[μ,:,:], shifts)
	end

	function substitute_U!(a::Abstractfield, b::Abstractfield)
		NX, NT = size(a)

        for it in 1:NT
            for ix in 1:NX
                for μ in 1:2
					a[μ,ix,it] = b[μ,ix,it]
				end
			end
		end

		return nothing
	end

    function add!(a::Abstractfield, b::Abstractfield, fac)
        NX, NT = size(a)

        for it in 1:NT
            for ix in 1:NX
                for μ in 1:2
					a[μ,ix,it] += fac * b[μ,ix,it]
				end
			end
		end

		return nothing
	end

	function random_gaugefield!(u::Gaugefield, rng)
		NX, NT = size(u)

		for it in 1:NT
			for ix in 1:NX
                for μ in 1:2
				    u[μ,ix,it] = (rand(rng) - 0.5) * 2 * 2π
                end
			end
		end

		recalc_Sg!(u)
		recalc_CV!(u)
		return nothing
	end

	function calc_Sg(u::Gaugefield)
		s = u.NV - plaquette_sum(u)
		return u.β * s
	end

	function recalc_Sg!(u::Gaugefield)
		u.Sg = calc_Sg(u)
		return nothing
	end

	function calc_CV(u::Gaugefield)
		NX, NT = size(u)
		q = 0.0

		for it in 1:NT; it_plu = mod1(it + 1, NT)
			for ix in 1:NX; ix_plu = mod1(ix + 1, NX)
				q += sin(u[1,ix,it] + u[2,ix_plu,it] - u[1,ix,it_plu] - u[2,ix,it])
			end
		end

		return q / 2π
	end

	function recalc_CV!(u::Gaugefield)
		u.CV = calc_CV(u)
		return nothing
	end

	function plaquette(u::Gaugefield, ix, it)
		NX, NT = size(u)
		ix_plu = mod1(ix + 1, NX)
		it_plu = mod1(it + 1, NT)
		return u[1,ix,it] + u[2,ix_plu,it] - u[1,ix,it_plu] - u[2,ix,it]
	end

	function plaquette_sum(u::Gaugefield)
		NX, NT = size(u)
		plaq = 0.0

		for it in 1:NT; it_plu = mod1(it + 1, NT)
			for ix in 1:NX; ix_plu = mod1(ix + 1, NX)
				plaq += cos(u[1,ix,it] + u[2,ix_plu,it] - u[1,ix,it_plu] - u[2,ix,it])
			end
		end

		return plaq
	end

	function wilson_loop(u::Gaugefield, LX, LT, ix, it; tempgauge=false)
		NX, NT = size(u)
		x_up_side = 0.0
		t_up_side = 0.0
		x_down_side = 0.0
		t_down_side = 0.0

		if tempgauge == true

			for i in 0:LX-1
				x_up_side += u[1,mod1(ix+i, NX),it]
				x_down_side -= u[1,mod1(ix+i, NX),mod1(it+LT, NT)]
			end

			return cis(x_up_side + x_down_side)
		elseif tempgauge == false

			for i in 0:LX-1
				x_up_side += u[1,mod1(ix+i, NX),it]
				x_down_side -= u[1,mod1(ix+i, NX),mod1(it+LT, NT)]
			end

			for i in 0:LT-1
				t_up_side += u[2,mod1(ix+LX, NX),mod1(it+i, NT)]
				t_down_side -= u[2,ix,mod1(it+i, NT)]
			end

			return cis(x_up_side + t_up_side + x_down_side + t_down_side)
		else
			error("tempgauge in Wilson loop calculation can only either be true or false")
		end
	end

	function staple(u::Gaugefield, μ, ix, it)
		NX, NT = size(u)
		it_min = mod1(it-1, NT)
		it_plu = mod1(it+1, NT)
		ix_min = mod1(ix-1, NX)
		ix_plu = mod1(ix+1, NX)

		if μ == 1
			l =  u[2,ix,it    ] + u[1,ix,it_plu] - u[2,ix_plu,it    ]
			r = -u[2,ix,it_min] + u[1,ix,it_min] + u[2,ix_plu,it_min]
		elseif μ == 2
			l =  u[1,ix    ,it] + u[2,ix_plu,it] - u[1,ix    ,it_plu]
			r = -u[1,ix_min,it] + u[2,ix_min,it] + u[1,ix_min,it_plu]
		end

		return cis(l) + cis(r)
	end

    include("liefields.jl")
end
