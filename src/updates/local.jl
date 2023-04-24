module Local
	using Random
	
	import ..System_parameters: Params
	import ..Gaugefields: Gaugefield,dqar,daction,add_Sg!,add_CV!
	import ..Metadynamics: Bias_potential,update_bias!,DeltaV
	import ..Verbose_print: Verbose_,println_verbose
	
	function metropolis!(field::Gaugefield, ix::Int64, it::Int64, μ::Int64, rng::Xoshiro, ϵ::Float64, multi_hit)
		accept = 0
		for hit = 1:multi_hit
			dU = randn(rng)*ϵ
			ΔSg = daction(field, ix, it, μ, dU)
			ΔCV = dqar(field, ix, it, μ, dU)
			if rand(rng) ≤ exp(-ΔSg)  
				@inbounds field[ix,it,μ] += dU
				add_Sg!(field, ΔSg)
				add_CV!(field, ΔCV)
				accept += 1
			end
		end
		return accept
	end
	
	function metropolis_meta!(field::Gaugefield, bias::Bias_potential, ix::Int64, it::Int64, μ::Int64, rng::Xoshiro, ϵ::Float64, multi_hit)
		accept = 0
		for hit = 1:multi_hit
			dU = randn(rng)*ϵ
			ΔCV = dqar(field, ix, it, μ, dU)

			oldCV = field.CV
			newCV = field.CV+ΔCV 
			
			ΔSg = daction(field, ix, it, μ, dU) 
			ΔV = DeltaV(bias, oldCV, newCV)
			if rand(rng) ≤ exp(-ΔSg-ΔV) 
				@inbounds field[ix,it,μ] += dU
				add_Sg!(field, ΔSg)
				add_CV!(field, ΔCV)
				accept += 1
			end
		end
		return accept
	end

	function sweep!(
		field::Gaugefield,
		rng::Xoshiro,
		ϵ::Float64,
		multi_hit,
		)
		NX,NT,_ = size(field)
		numaccepts = 0
		for eo in [0,1]
			for ix = (eo+1):2:NX
				for it = 1:NT
					accept = metropolis!(field, ix, it, 2, rng, ϵ, multi_hit)
					numaccepts += accept
				end
			end
		end
		for eo in [0,1]
			for it = (eo+1):2:NT
				for ix = 1:NX
					accept = metropolis!(field, ix, it, 1, rng, ϵ, multi_hit)
					numaccepts += accept
				end
			end
		end
		return numaccepts
	end

	function sweep_meta!(
		field::Gaugefield,
		bias::Bias_potential,
		rng::Xoshiro,
		ϵ::Float64,
		multi_hit,
		)
		NX,NT,_ = size(field)
		numaccepts = 0
		for eo in [0,1]
			for ix = (eo+1):2:NX
				for it = 1:NT
					accept = metropolis_meta!(field, bias, ix, it, 2, rng, ϵ, multi_hit)
					numaccepts += accept
				end
			end
		end
		for eo in [0,1]
			for it = (eo+1):2:NT
				for ix = 1:NX
					accept = metropolis_meta!(field, bias, ix, it, 1, rng, ϵ, multi_hit)
					numaccepts += accept
				end
			end
		end
		return numaccepts
	end

	function adjusted_ϵ(ϵ::Float64, numaccepts, metro_norm, metro_target_acc)
		return 	ϵ + (numaccepts*metro_norm - metro_target_acc) * 0.2
	end

end