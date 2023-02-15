module MC
	using Random
	using Printf
	#using Base.Threads:@threads,nthreads,threadid
	
	import ..System_parameters:Params
	import ..Gaugefields:Gaugefield,recalc_S!,recalc_Q!,dqar,daction,plaquette,swap!
	import ..Metadynamics:Bias_potential,update_bias!,penalty_potential
	import ..Measurements:Measurement_set,measurements
	import ..Verbose_print:Verbose_,println_verbose

	export sweep!,sweep_meta!,try_swap!
	
	function metropolis!(field::Gaugefield,ix::Int64,it::Int64,μ::Int64,rng::Xoshiro,ϵ::Float64)
		dU = randn(rng)*ϵ
		ΔS = daction(field,ix,it,μ,dU)
		ΔQ = dqar(field,ix,it,μ,dU)
		accept = rand(rng) ≤ exp(-ΔS)
		if accept  
			@inbounds field.g[ix,it,μ] += dU
			field.Q += ΔQ 
			field.S += ΔS
		end
		return accept
	end
	
	function metropolis_meta!(field::Gaugefield,bias::Bias_potential,ix::Int64,it::Int64,μ::Int64,rng::Xoshiro,ϵ::Float64,static::Bool)
		@inline index(q_val,q_min,dq) = round(Int,(q_val-q_min)/dq+0.5,RoundNearestTiesAway)
		dU = randn(rng)*ϵ
		ΔQ = dqar(field,ix,it,μ,dU)
		
		old_ind = index(field.Q,bias.Qmin,bias.δq)
		prop_ind = index(field.Q+ΔQ,bias.Qmin,bias.δq)
		
		ΔS = daction(field,ix,it,μ,dU) 
		ΔV = bias[prop_ind]-bias[old_ind]
		ΔV_pen = penalty_potential(field.Q+ΔQ,bias.Qmin_thr,bias.Qmax_thr,bias.k)-penalty_potential(field.Q,bias.Qmin_thr,bias.Qmax_thr,bias.k)
		accept = rand(rng) ≤ exp(-ΔS-ΔV-ΔV_pen) 
		if accept 
			@inbounds field.g[ix,it,μ] += dU
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
			for ix = (eo+1):2:field.NX
				for it = 1:field.NT
					accept = metropolis!(field,ix,it,2,rng,ϵ)
					numaccepts += ifelse(accept,1,0)
				end
			end
		end
		for eo in [0,1]
			for it = (eo+1):2:field.NT
				for ix = 1:field.NX
					accept = metropolis!(field,ix,it,1,rng,ϵ)
					numaccepts += ifelse(accept,1,0)
				end
			end
		end
		return numaccepts
	end

	function sweep_meta!(field::Gaugefield,bias::Bias_potential,rng::Xoshiro,ϵ::Float64,static::Bool)
		numaccepts = 0
		for eo in [0,1]
			for ix = (eo+1):2:field.NX
				for it = 1:field.NT
					accept = metropolis_meta!(field,bias,ix,it,2,rng,ϵ,static)
					numaccepts += ifelse(accept,1,0)
				end
			end
		end
		for eo in [0,1]
			for it = (eo+1):2:field.NT
				for ix = 1:field.NX
					accept = metropolis_meta!(field,bias,ix,it,1,rng,ϵ,static)
					numaccepts += ifelse(accept,1,0)
				end
			end
		end
		return sum(numaccepts)
	end

	function try_swap!(field_main::Gaugefield,field_meta::Gaugefield,bias::Bias_potential,rng::Xoshiro,static::Bool)
		@inline index(q_val,q_min,dq) = round(Int,(q_val-q_min)/dq+0.5,RoundNearestTiesAway)
		main_ind = index(field_main.Q,bias.Qmin,bias.δq)
		meta_ind = index(field_meta.Q,bias.Qmin,bias.δq)

		ΔQ = field_meta.Q - field_main.Q
		ΔS = field_meta.S - field_main.S
		ΔV = bias[meta_ind] - bias[main_ind]
		ΔV_pen = penalty_potential(field_meta.Q,bias.Qmin_thr,bias.Qmax_thr,bias.k) - penalty_potential(field_main.Q,bias.Qmin_thr,bias.Qmax_thr,bias.k)
		accept_swap = rand(rng) ≤ exp(ΔV+ΔV_pen)
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

end