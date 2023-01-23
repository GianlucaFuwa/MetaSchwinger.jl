module MC
	using Random
	using Printf
	
	import ..System_parameters:Params
	import ..Gaugefields:Gaugefield,recalc_S!,recalc_Q!,dqar,daction,plaquette,swap!
	import ..Metadynamics:Bias_potential,update_bias!,penalty_potential
	import ..Measurements:Measurement_set,measurements
	import ..Verbose_print:Verbose_,println_verbose
	
	function metropolis!(field::Gaugefield,nx::Int64,nt::Int64,d::Int64,rng::Xoshiro,ϵ::Float64)
		dU = (rand(rng)-0.5)*2*ϵ
		ΔS = daction(field,nx,nt,d,dU)
		ΔQ = dqar(field,nx,nt,d,dU)
		accept = rand(rng) ≤ exp(-ΔS)
		if accept  
			field.g[nx,nt,d] += dU
			field.Q += ΔQ 
			field.S += ΔS
		end
		return accept
	end
	
	function metropolis_meta!(field::Gaugefield,bias::Bias_potential,nx::Int64,nt::Int64,d::Int64,rng::Xoshiro,ϵ::Float64,static::Bool)
		@inline index(q_val,q_min,dq) = round(Int,(q_val-q_min)/dq+0.5,RoundNearestTiesAway)
		dU = (rand(rng)-0.5)*2*ϵ
		ΔQ = dqar(field,nx,nt,d,dU)
		
		old_ind = index(field.Q,bias.Qmin,bias.δq)
		prop_ind = index(field.Q+ΔQ,bias.Qmin,bias.δq)
		
		ΔS = daction(field,nx,nt,d,dU) 
		ΔV = bias.values[prop_ind]-bias.values[old_ind]
		ΔV_pen = penalty_potential(field.Q+ΔQ,bias.Qmin_thr,bias.Qmax_thr,bias.k)-penalty_potential(field.Q,bias.Qmin_thr,bias.Qmax_thr,bias.k)
		accept = rand(rng) ≤ exp(-ΔS-ΔV-ΔV_pen) 
		if accept 
			field.g[nx,nt,d] += dU
			field.Q += ΔQ 
			field.S += ΔS
			if ~static 
				update_bias!(bias,field.Q)
			end
		end
		return accept
	end

	function sweep!(field::Gaugefield,rng::Xoshiro,ϵ::Float64)
		numaccepts = 0
		for eo in [0,1]
			for nx = (eo+1):2:field.Nx
				for nt = 1:field.Nt
					accept = metropolis!(field,nx,nt,2,rng,ϵ)
					numaccepts += ifelse(accept,1,0)
				end
			end
		end
		for eo in [0,1]
			for nt = (eo+1):2:field.Nt
				for nx = 1:field.Nx
					accept = metropolis!(field,nx,nt,1,rng,ϵ)
					numaccepts += ifelse(accept,1,0)
				end
			end
		end
		return numaccepts
	end

	function sweep_meta!(field::Gaugefield,bias::Bias_potential,rng::Xoshiro,ϵ::Float64,static::Bool)
		numaccepts = 0
		for eo in [0,1]
			for nx = (eo+1):2:field.Nx
				for nt = 1:field.Nt
					accept = metropolis_meta!(field,bias,nx,nt,2,rng,ϵ,static)
					numaccepts += ifelse(accept,1,0)
				end
			end
		end
		for eo in [0,1]
			for nt = (eo+1):2:field.Nt
				for nx = 1:field.Nx
					accept = metropolis_meta!(field,bias,nx,nt,1,rng,ϵ,static)
					numaccepts += ifelse(accept,1,0)
				end
			end
		end
		return numaccepts
	end

	function try_swap!(field_main::Gaugefield,field_meta::Gaugefield,bias::Bias_potential,rng::Xoshiro,static::Bool,verbose::Verbose_)
		@inline index(q_val,q_min,dq) = round(Int,(q_val-q_min)/dq+0.5,RoundNearestTiesAway)
		old_ind = index(field_main.Q,bias.Qmin,bias.δq)
		prop_ind = index(field_meta.Q,bias.Qmin,bias.δq)

		ΔQ = field_meta.Q - field_main.Q
		ΔS = field_meta.S - field_main.S
		ΔV = bias.values[prop_ind] - bias.values[old_ind]
		ΔV_pen = penalty_potential(field_meta.Q,bias.Qmin_thr,bias.Qmax_thr,bias.k) - penalty_potential(field_main.Q,bias.Qmin_thr,bias.Qmax_thr,bias.k)
		accept_swap = rand(rng) ≤ exp(-ΔV-ΔV_pen)
		if accept_swap
			swap!(field_main,field_meta)
			field_meta.Q -= ΔQ
			field_main.Q += ΔQ
			field_meta.S -= ΔS
			field_main.S += ΔS
			if ~static
				update_bias!(bias,field_meta.Q)
			end
		end
		return accept_swap
	end 


	function instanton_update!(field::Gaugefield,bias::Bias_potential,rng::Xoshiro,insta::Array{Float64,3})
		@inline index(q_val,q_min,dq) = round(Int,(q_val-q_min)/dq+0.5,RoundNearestTiesAway)
		Q = deepcopy(field.Q)
		S = deepcopy(field.S)
		plus_minus = sign(rand(rng)-0.5)
		field.g .+= plus_minus * insta

		recalc_Q!(field)
		recalc_S!(field)
		ΔQ = field.Q-Q
		ΔS = field.S-S

		old_ind = index(Q,bias.Qmin,bias.δq)
		prop_ind = index(Q+ΔQ,bias.Qmin,bias.δq)
		
		ΔV = bias.values[prop_ind]-bias.values[old_ind]
		ΔV_pen = penalty_potential(Q+ΔQ,bias.Qmin_thr,bias.Qmax_thr,bias.k)-penalty_potential(Q,bias.Qmin_thr,bias.Qmax_thr,bias.k)
		accept = rand(rng) ≤ exp(-ΔS-ΔV-ΔV_pen)
		if accept   
			update_bias!(bias,field.Q)
		else
			field.g .-= plus_minus * insta
		end
		return accept
	end

end