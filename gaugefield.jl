module Gaugefields
	import ..System_parameters:Params
	
	mutable struct Gaugefield
		NX::Int64
		NT::Int64
		NV::Int64
		β::Float64
		g::Array{Float64,3}		
		Q::Float64
		S::Float64
		
		function Gaugefield(p::Params)
			NX = p.N[1]
			NT = p.N[2]
			NV = NX*NT
			g = zeros(Float64,NX,NT,2)
			Q = 0.
			S = 0.
			return new(NX,NT,NV,p.β,g,Q,S)
		end
	end
	
	function recalc_Q!(x::Gaugefield)
		NX = x.NX
		NT = x.NT
		q = 0.
		for it=1:NT
			for ix=1:NX
			q += sin(plaquette(x.g,ix,it))
			end
		end
		x.Q = q / 2pi
		return nothing
	end
	
	function recalc_S!(x::Gaugefield)
		NX = x.NX
		NT = x.NT
		s = 0.
		for it=1:NT
			for ix=1:NX
				s += 1-cos(plaquette(x.g,ix,it))
			end
		end
		x.S = x.β*s
		return nothing
	end
	
	function plaquette(g::Array{Float64,3},ix,it)
		NX,NT, = size(g)
		return g[ix,it,1]+g[mod1(ix+1,NX),it,2]-g[ix,mod1(it+1,NT),1]-g[ix,it,2]
	end
	
	function Wilson_loop(g::Array{Float64,3},x,t,ix,it)
		NX,NT, = size(g)
		x_up_side = 0.
		t_up_side = 0.
		x_down_side = 0.
		t_down_side = 0.
		for i=0:x-1
			x_up_side += g[mod1(ix+i,NX),it,1] 
			x_down_side -= g[mod1(ix+i,NX),mod1(it+t,NT),1]
		end
		for i=0:t-1
			t_up_side += g[mod1(ix+x,NX),mod1(it+i,NT),2]
			t_down_side -= g[ix,mod1(it+i,NT),2]
		end
		return x_up_side + t_up_side + x_down_side + t_down_side
	end
	
	function staple(links::Array{Float64,3},ix::Int64,it::Int64,μ::Int64)
		NX = size(links,1)
		@inline →(n) = mod1(n+1,NX) 
		@inline ←(n) = mod1(n-1,NX) 
		if μ == 1
			l = links[ix,it,2]     + links[ix,→(it),1] - links[→(ix),it,2]
			r = -links[ix,←(it),2] + links[ix,←(it),1] + links[→(ix),←(it),2]
		elseif μ == 2
			l = links[ix,it,1]     + links[→(ix),it,2] - links[ix,→(it),1]
			r = -links[←(ix),it,1] + links[←(ix),it,2] + links[←(ix),→(it),1]
		end
		return exp(l*im)+exp(r*im)
	end
	
	function daction(x::Gaugefield,ix::Int64,it::Int64,μ::Int64,dU::Float64)
		link_old = x.g[ix,it,μ]
		A = staple(x.g,ix,it,μ)
		return -x.β*real( (exp((link_old+dU)*im) - exp(link_old*im)) * A' ) 
	end
	
	function dqar(x::Gaugefield,ix::Int64,it::Int64,μ::Int64,dU::Float64)
		links = x.g
		@inline →(n) = mod1(n+1,x.NX) 
		@inline ←(n) = mod1(n-1,x.NX)   

		if μ == 1
			a = links[ix,it,1]+dU + links[→(ix),it,2]    - links[ix,→(it),1] - links[ix,it,2]
			b = links[ix,←(it),1] + links[→(ix),←(it),2] - links[ix,it,1]-dU - links[ix,←(it),2]
			c = links[ix,it,1]    + links[→(ix),it,2]    - links[ix,→(it),1] - links[ix,it,2]
			d = links[ix,←(it),1] + links[→(ix),←(it),2] - links[ix,it,1]    - links[ix,←(it),2]
		elseif μ == 2
			a = links[←(ix),it,1] + links[ix,it,2]+dU - links[←(ix),→(it),1] - links[←(ix),it,2]
			b = links[ix,it,1]    + links[→(ix),it,2] - links[ix,→(it),1]    - links[ix,it,2]-dU
			c = links[←(ix),it,1] + links[ix,it,2]    - links[←(ix),→(it),1] - links[←(ix),it,2]
			d = links[ix,it,1]    + links[→(ix),it,2] - links[ix,→(it),1]    - links[ix,it,2]
		end

		return (sin(a) + sin(b) - sin(c) - sin(d)) / 2pi
	end
	
	function swap!(a::Gaugefield,b::Gaugefield)
		tmp = deepcopy(a.g)
		a.g = deepcopy(b.g)
		b.g = tmp
		return nothing
	end
	
end
		
		