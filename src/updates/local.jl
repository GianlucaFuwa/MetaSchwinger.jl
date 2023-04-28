module Local
	using Random
	
	import ..Gaugefields: add_CV!, add_Sg!, daction, dqar, Gaugefield
	import ..Metadynamics: Bias_potential, DeltaV, update_bias!
	import ..System_parameters: Params
	import ..Verbose_print: println_verbose, Verbose_
	
	function metropolis!(field::Gaugefield, ix, it, μ, rng, ϵ, multi_hit)
		accept = 0

		for hit in 1:multi_hit
			dU = randn(rng) * ϵ
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
	
	function metropolis_meta!(field::Gaugefield, bias::Bias_potential, ix, it, μ, rng, ϵ, multi_hit)
		accept = 0

		for hit in 1:multi_hit
			dU = randn(rng) * ϵ
			ΔCV = dqar(field, ix, it, μ, dU)

			oldCV = field.CV
			newCV = field.CV + ΔCV 
			
			ΔSg = daction(field, ix, it, μ, dU) 
			ΔV = DeltaV(bias, oldCV, newCV)

			if rand(rng) ≤ exp(-ΔSg - ΔV) 
				@inbounds field[ix,it,μ] += dU
				add_Sg!(field, ΔSg)
				add_CV!(field, ΔCV)
				accept += 1
			end
		end

		return accept
	end

	function sweep!(field::Gaugefield, rng, ϵ, multi_hit)
		NX, NT, _ = size(field)
		numaccepts = 0

		for eo in [0,1]
			for ix in (eo+1):2:NX
				for it in 1:NT
					accept = metropolis!(field, ix, it, 2, rng, ϵ, multi_hit)
					numaccepts += accept
				end
			end
		end

		for eo in [0,1]
			for it in (eo+1):2:NT
				for ix in 1:NX
					accept = metropolis!(field, ix, it, 1, rng, ϵ, multi_hit)
					numaccepts += accept
				end
			end
		end

		return numaccepts
	end

	function sweep_meta!(field::Gaugefield, bias::Bias_potential, rng, ϵ, multi_hit)
		NX, NT, _ = size(field)
		numaccepts = 0

		for eo in [0,1]
			for ix in (eo+1):2:NX
				for it in 1:NT
					accept = metropolis_meta!(field, bias, ix, it, 2, rng, ϵ, multi_hit)
					numaccepts += accept
				end
			end
		end

		for eo in [0,1]
			for it in (eo+1):2:NT
				for ix in 1:NX
					accept = metropolis_meta!(field, bias, ix, it, 1, rng, ϵ, multi_hit)
					numaccepts += accept
				end
			end
		end

		return numaccepts
	end

	function adjusted_ϵ(ϵ, numaccepts, metro_norm, metro_target_acc)
		return 	ϵ + (numaccepts * metro_norm - metro_target_acc) * 0.2
	end

end