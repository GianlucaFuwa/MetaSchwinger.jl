module Gaugefields
	import ..System_parameters:Params
	export daction
	export dqar
	export recalc_S!
	export recalc_Q!
	export plaquette
	export staple
	export instanton
	
	mutable struct Gaugefield
		Nx::Int64
		Nt::Int64
		Nd::Int64
		β::Float64
		g::Array{Float64,3}		
		Q::Float64
		S::Float64
		
		function Gaugefield(p::Params)
			Nx = p.N[1]
			Nt = p.N[2]
			Nd = 2
			g = zeros(Float64,Nx,Nt,Nd)
			Q = 0.
			S = 0.
			return new(Nx,Nt,Nd,p.β,g,Q,S)
		end
	end
	
	function recalc_Q!(x::Gaugefield)
		q = 0. 
		for nx=1:x.Nx, nt=1:x.Nt
			q += sin(plaquette(x.g,nx,nt))
		end
		x.Q = q / 2pi
		return nothing
	end
	
	function recalc_S!(x::Gaugefield)
		s = 0.
		for nx=1:x.Nx, nt=1:x.Nt
			s += 1-cos(plaquette(x.g,nx,nt))
		end
		x.S = x.β*s
		return nothing
	end
	
	function plaquette(g::Array{Float64,3},nx,nt)
		Nx,Nt, = size(g)
		return g[nx,nt,1]+g[mod1(nx+1,Nx),nt,2]-g[nx,mod1(nt+1,Nt),1]-g[nx,nt,2]
	end
	
	function Wilson_loop(g::Array{Float64,3},x,t,nx,nt)
		Nx,Nt, = size(g)
		x_up_side = 0.
		t_up_side = 0.
		x_down_side = 0.
		t_down_side = 0.
		for i=0:x-1
			x_up_side += g[mod1(nx+i,Nx),nt,1] 
			x_down_side -= g[mod1(nx+i,Nx),mod1(nt+t,Nt),1]
		end
		for i=0:t-1
			t_up_side += g[mod1(nx+x,Nx),mod1(nt+i,Nt),2]
			t_down_side -= g[nx,mod1(nt+i,Nt),2]
		end
		return x_up_side + t_up_side + x_down_side + t_down_side
	end
	
	function staple(links::Array{Float64,3},nx::Int64,nt::Int64,μ::Int64)
		Nx = size(links,1)
		@inline →(n) = mod1(n+1,Nx) 
		@inline ←(n) = mod1(n-1,Nx) 
		if μ == 1
			l = links[nx,nt,2]     + links[nx,→(nt),1] - links[→(nx),nt,2]
			r = -links[nx,←(nt),2] + links[nx,←(nt),1] + links[→(nx),←(nt),2]
		elseif μ == 2
			l = links[nx,nt,1]     + links[→(nx),nt,2] - links[nx,→(nt),1]
			r = -links[←(nx),nt,1] + links[←(nx),nt,2] + links[←(nx),→(nt),1]
		end
		return exp(l*im)+exp(r*im)
	end
	
	function daction(x::Gaugefield,nx::Int64,nt::Int64,μ::Int64,dU::Float64)
		link_old = x.g[nx,nt,μ]
		A = staple(x.g,nx,nt,μ)
		return -x.β*real( (exp((link_old+dU)*im) - exp(link_old*im))*conj(A) ) 
	end
	
	function dqar(x::Gaugefield,nx::Int64,nt::Int64,μ::Int64,dU::Float64)
		links = x.g
		@inline →(n) = mod1(n+1,x.Nx) 
		@inline ←(n) = mod1(n-1,x.Nx)   

		if μ == 1
			a = links[nx,nt,1]+dU + links[→(nx),nt,2]    - links[nx,→(nt),1] - links[nx,nt,2]
			b = links[nx,←(nt),1] + links[→(nx),←(nt),2] - links[nx,nt,1]-dU - links[nx,←(nt),2]
			c = links[nx,nt,1]    + links[→(nx),nt,2]    - links[nx,→(nt),1] - links[nx,nt,2]
			d = links[nx,←(nt),1] + links[→(nx),←(nt),2] - links[nx,nt,1]    - links[nx,←(nt),2]
		elseif μ == 2
			a = links[←(nx),nt,1] + links[nx,nt,2]+dU - links[←(nx),→(nt),1] - links[←(nx),nt,2]
			b = links[nx,nt,1]    + links[→(nx),nt,2] - links[nx,→(nt),1]    - links[nx,nt,2]-dU
			c = links[←(nx),nt,1] + links[nx,nt,2]    - links[←(nx),→(nt),1] - links[←(nx),nt,2]
			d = links[nx,nt,1]    + links[→(nx),nt,2] - links[nx,→(nt),1]    - links[nx,nt,2]
		end

		return (sin(a) + sin(b) - sin(c) - sin(d)) / 2pi
	end

	function instanton(x::Gaugefield,n::Int64)
		n_insta = zeros(size(x.g))
		for nx = 1:x.Nx
			n_insta[nx,x.Nt,2] += n*nx*2pi/x.Nx
		end
		for nt = 1:x.Nt
			n_insta[:,nt,1] .+= -ones(x.Nx).*n*nt*2pi/x.Nx/x.Nt 
		end
		return n_insta
	end

	function swap!(a::Gaugefield,b::Gaugefield)
		@assert size(a.g) == size(b.g)
		tmp = deepcopy(a.g)
		a.g = deepcopy(b.g)
		b.g = tmp
		return nothing
	end
	
end
		
		