module Metadynamics
	using ..System_parameters:Params
	
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
		fp::IOStream
		
		function Bias_potential(p::Params,instance::Int64)
			Qmin = p.Qmax[1]
			Qmax = p.Qmax[2]
			Qmin_thr = p.Qthr[1]
			Qmax_thr = p.Qthr[2]
			δq = p.δq
			w = p.w
			k = p.k
			values = potential_from_file(p,p.usebiases[instance])
			q_vals = range(p.Qmax[1],p.Qmax[2],Int((p.Qmax[2]-p.Qmax[1])/p.δq)+1)
			fp = open(p.biasfiles[instance],"w")
			return new(Qmin,Qmax,Qmin_thr,Qmax_thr,δq,w,k,values,q_vals,fp)
		end	
	end

	function potential_from_file(p::Params,usebias::Union{Nothing,IOStream})
		if usebias === nothing
			return zeros(round(Int,(p.Qmax[2]-p.Qmax[1])/p.δq,RoundNearestTiesAway)+1)
		else
			values = readdlm(usebias,Float64)
			@assert length(values[:,2]) == round(Int,((p.Qmax[2]-p.Qmax[1])/p.δq),RoundNearestTiesAway)+1 "Potential doesn't fit its parameters"
			return values[:,2]
		end
	end

	function Base.flush(x::Bias_potential)
		if x.fp !== nothing
		    flush(x.fp)
		end
    	end

	function Base.setindex!(x::Bias_potential,v,i::Int)
		x.values[min(length(x.values),max(i,1))] = v
		return nothing
	end

	function Base.getindex(x::Bias_potential,i::Int)
		return x.values[min(length(x.values),max(i,1))]
	end

	function grid_ind(q::Float64,q_min::Int64,δq::Float64)
		grid_index = (q-q_min)/δq + 0.5
		return round(Int,grid_index,RoundNearestTiesAway)
	end

	function update_bias!(x::Bias_potential,q::Float64)
		@inline index(q_val,q_min,dq) = round(Int,(q_val-q_min)/dq+0.5,RoundNearestTiesAway)
		grid_index = index(q,x.Qmin,x.δq)
		grid_q = x.Qmin + grid_index*x.δq

		x[grid_index] += x.w*(1-(q-grid_q)/x.δq)
		x[grid_index+1] += x.w*(q-grid_q)/x.δq
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

end


