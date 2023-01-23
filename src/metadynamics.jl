module Metadynamics
	using Printf
	using ..System_parameters:Params
	export penalty_potential
	export update_bias!
	export Bias_potential
	export make_symmetric!
	
	mutable struct Bias_potential
		Qmin::Float64
		Qmax::Float64
		Qmin_thr::Float64
		Qmax_thr::Float64
		δq::Float64
		w::Float64
		k::Int64
		
		values::Array{Float64,1}
		q_vals::Array{Float64,1}
		
		function Bias_potential(p::Params)
			Qmin = p.Qmax[1]
			Qmax = p.Qmax[2]
			Qmin_thr = p.Qthr[1]
			Qmax_thr = p.Qthr[2]
			δq = p.δq
			w = p.w
			k = p.k
			values = p.usebias
			q_vals = range(p.Qmax[1],p.Qmax[2],Int((p.Qmax[2]-p.Qmax[1])/p.δq)+1)
			return new(Qmin,Qmax,Qmin_thr,Qmax_thr,δq,w,k,values,q_vals)
		end	
	end
	
	function grid_ind(q::Float64,q_min::Int64,δq::Float64)
		grid_index = (q-q_min)/δq + 0.5
		return round(Int,grid_index,RoundNearestTiesAway)
	end

	function update_bias!(x::Bias_potential,q::Float64)
		@inline index(q_val,q_min,dq) = round(Int,(q_val-q_min)/dq+0.5,RoundNearestTiesAway)
		grid_index = index(q,x.Qmin,x.δq)
		grid_q = x.q_vals[grid_index]

		x.values[grid_index] += x.w*(1-(q-grid_q)/x.δq)
		x.values[grid_index+1] += x.w*(q-grid_q)/x.δq
		return nothing
	end

	function penalty_potential(q::Float64,Qmin_thr::Float64,Qmax_thr::Float64,k)
		if q < Qmin_thr || q > Qmax_thr
			p_pot = k*min((q-Qmin_thr)^2,(q-Qmax_thr)^2)
		else 
			p_pot = 0
		end
		return p_pot
	end

	function make_symmetric!(b::Bias_potential)
		mid = ceil(Int,length(b.values)/2)
		for idx = 1:mid-1
			idx_mir = -idx+2*mid
			b.values[idx] = (b.values[idx]+b.values[idx_mir])/2
			b.values[idx_mir] = b.values[idx]
		end
		return nothing
	end
end


