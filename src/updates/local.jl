module Local
	using Random
	
	import ..System_parameters: Params
	import ..Gaugefields: Gaugefield,dqar,daction,add_Sg!,add_CV!
	import ..Metadynamics: Bias_potential,update_bias!,DeltaV
	import ..Verbose_print: Verbose_,println_verbose
	
	function metropolis!(field::Gaugefield, ix::Int64, it::Int64, μ::Int64, rng::Xoshiro, ϵ::Float64)
		dU = randn(rng)*ϵ
		ΔSg = daction(field, ix, it, μ, dU)
		ΔCV = dqar(field, ix, it, μ, dU)
		accept = rand(rng) ≤ exp(-ΔSg)
		if accept  
			@inbounds field[ix,it,μ] += dU
			add_Sg!(field, ΔSg)
			add_CV!(field, ΔCV)
		end
		return accept
	end
	
	function metropolis_meta!(field::Gaugefield, bias::Bias_potential, ix::Int64, it::Int64, μ::Int64, rng::Xoshiro, ϵ::Float64)
		dU = randn(rng)*ϵ
		ΔCV = dqar(field, ix, it, μ, dU)

		oldCV = field.CV
		newCV = field.CV+ΔCV 
		
		ΔSg = daction(field, ix, it, μ, dU) 
		ΔV = DeltaV(bias, oldCV, newCV)
		accept = rand(rng) ≤ exp(-ΔSg-ΔV) 
		if accept 
			@inbounds field[ix,it,μ] += dU
			add_Sg!(field, ΔSg)
			add_CV!(field, ΔCV)
		end
		return accept
	end

	function sweep!(field::Gaugefield, rng::Xoshiro, ϵ::Float64)
		NX,NT,_ = size(field)
		numaccepts = 0
		for eo in [0,1]
			for ix = (eo+1):2:NX
				for it = 1:NT
					accept = metropolis!(field, ix, it, 2, rng, ϵ)
					numaccepts += ifelse(accept,1,0)
				end
			end
		end
		for eo in [0,1]
			for it = (eo+1):2:NT
				for ix = 1:NX
					accept = metropolis!(field, ix, it, 1, rng, ϵ)
					numaccepts += ifelse(accept,1,0)
				end
			end
		end
		return numaccepts
	end

	function sweep_meta!(field::Gaugefield, bias::Bias_potential, rng::Xoshiro, ϵ::Float64)
		NX,NT,_ = size(field)
		numaccepts = 0
		for eo in [0,1]
			for ix = (eo+1):2:NX
				for it = 1:NT
					accept = metropolis_meta!(field, bias, ix, it, 2, rng, ϵ)
					numaccepts += ifelse(accept, 1, 0)
				end
			end
		end
		for eo in [0,1]
			for it = (eo+1):2:NT
				for ix = 1:NX
					accept = metropolis_meta!(field, bias, ix, it, 1, rng, ϵ)
					numaccepts += ifelse(accept, 1, 0)
				end
			end
		end
		return sum(numaccepts)
	end

end