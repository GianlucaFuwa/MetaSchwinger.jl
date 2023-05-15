module Gaugefields

	import ..System_parameters: Params
	
	struct Gaugefield
		NX::Int64
		NT::Int64
		NV::Int64
		β::Float64
		U::Array{Float64,3}		
		# RefValue lets us change Sg and CV eventhough Gaugefield and Floats are static
		Sg::Base.RefValue{Float64}
		CV::Base.RefValue{Float64}
		
		function Gaugefield(p::Params)
			NX = p.N[1]
			NT = p.N[2]
			NV = NX * NT
			β = p.β
			U = zeros(Float64, NX, NT, 2)
			Sg = Base.RefValue{Float64}(0.0)
			CV = Base.RefValue{Float64}(0.0)
			return new(NX, NT, NV, β, U, Sg, CV)
		end
	end

	function Base.setindex!(g::Gaugefield, v, ix, it, μ)
		@inbounds g.U[ix,it,μ] = v
	end

	@inline function Base.getindex(g::Gaugefield, ix, it, μ)
		@inbounds return g.U[ix,it,μ]
	end

	function Base.size(g::Gaugefield)
		return g.NX, g.NT, 2
	end
	# Need to overload get- and setproperty! to be able to change Sg and CV
	function Base.getproperty(g::Gaugefield, p::Symbol)
		if p === :Sg 
			return getfield(g, p)[]
		elseif p === :CV
			return getfield(g, p)[]
		else 
			return getfield(g, p)
		end
	end

	function Base.setproperty!(g::Gaugefield, p::Symbol, val)
		if p === :Sg 
			getfield(g, p)[] = val
		elseif p === :CV
			getfield(g, p)[] = val
		else 
			setfield!(g, p, val)
		end

		return nothing
	end	
	
	function Base.circshift(g::Gaugefield, shifts::NTuple{2,Int}, μ::Int64)
		return circshift(g[:,:,μ], shifts)
	end

	function substitute_U!(g::Gaugefield, U)
		NX, NT, d = size(g)

		for μ in 1:d
			for it in 1:NT
				for ix in 1:NX
					g[ix,it,μ] = U[ix,it,μ]
				end
			end
		end

		return nothing
	end

	function calc_Sg(g::Gaugefield)
		NX, NT, _ = size(g)
		s = 0.0

		for it in 1:NT
			for ix in 1:NX
				s += 1 - cos(plaquette(g, ix, it))
			end
		end

		return g.β * s
	end
	
	function recalc_Sg!(g::Gaugefield)
		g.Sg = calc_Sg(g)
		return nothing
	end

	function recalc_CV!(g::Gaugefield)
		NX, NT, _ = size(g)
		q = 0.0
		 
		for it in 1:NT
			for ix in 1:NX
				q += sin(plaquette(g, ix, it))
			end
		end

		g.CV = q / 2π
		return nothing
	end
	
	function plaquette(g::Gaugefield, ix, it)
		NX, NT, _ = size(g)
		return g[ix,it,1] + g[mod1(ix+1,NX),it,2] - g[ix,mod1(it+1,NT),1] - g[ix,it,2]
	end

	function plaquette_sum(g::Gaugefield)
		NX, NT ,_ = size(g)
		plaq_real = 0.0
		plaq_imag = 0.0

		for it in 1:NT
			for ix in 1:NX
				plaq_real += cos(plaquette(g, ix, it))
				plaq_imag += sin(plaquette(g, ix, it))
			end
		end

		return plaq_real, plaq_imag
	end
	
	function wilson_loop(g::Gaugefield, LX, LT, ix, it; tempgauge = false)
		NX, NT, _ = size(g)
		x_up_side = 0.0
		t_up_side = 0.0
		x_down_side = 0.0
		t_down_side = 0.0

		if tempgauge == true

			for i in 0:LX-1
				x_up_side += g[mod1(ix+i, NX),it,1] 
				x_down_side -= g[mod1(ix+i, NX),mod1(it+LT, NT),1]
			end
			
			return exp(x_up_side + x_down_side)
		elseif tempgauge == false

			for i in 0:LX-1
				x_up_side += g[mod1(ix+i, NX),it,1] 
				x_down_side -= g[mod1(ix+i, NX),mod1(it+LT, NT),1]
			end

			for i in 0:LT-1
				t_up_side += g[mod1(ix+LX, NX),mod1(it+i, NT),2]
				t_down_side -= g[ix,mod1(it+i, NT),2]
			end

			return x_up_side + t_up_side + x_down_side + t_down_side
		else
			error("tempgauge in Wilson loop calculation can only either be true or false")
		end
	end
	
	function staple(g::Gaugefield, ix, it, μ)
		NX, NT, _ = size(g)
		it_min = mod1(it-1, NT)
		it_plu = mod1(it+1, NT)
		ix_min = mod1(ix-1, NX)
		ix_plu = mod1(ix+1, NX)

		if μ == 1
			l =  g[ix,it    ,2] + g[ix,it_plu,1] - g[ix_plu,it    ,2]
			r = -g[ix,it_min,2] + g[ix,it_min,1] + g[ix_plu,it_min,2]
		elseif μ == 2
			l =  g[ix    ,it,1] + g[ix_plu,it,2] - g[ix    ,it_plu,1]
			r = -g[ix_min,it,1] + g[ix_min,it,2] + g[ix_min,it_plu,1]
		end

		return cis(l) + cis(r)
	end
	
end
		
		
