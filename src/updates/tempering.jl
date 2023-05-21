module Tempering
    
    import ..Gaugefields: Gaugefield
    import ..Metadynamics: BiasPotential, update_bias!

    function tempering_swap!(
		field1::Gaugefield,
		field2::Gaugefield,
		bias1::BiasPotential,
		bias2::BiasPotential,
		rng;
		actually_swap::Bool = true, # "Fake swap" for when we want to do parallel building
	)
		CV1 = field1.CV
		CV2 = field2.CV

		ΔCV = CV1 - CV2
		ΔSg = field1.Sg - field2.Sg
		ΔV1 = bias1(CV2) - bias1(CV1)
		ΔV2 = bias2(CV1) - bias2(CV2)
		accept_swap = rand(rng) ≤ exp(- ΔV1 - ΔV2)

		if accept_swap && actually_swap
			swap!(field1, field2)
			field2.CV += ΔCV
			field1.CV -= ΔCV
			field2.Sg += ΔSg
			field1.Sg -= ΔSg
			
			if !bias1.is_static
				update_bias!(bias1, field1.CV)
			end
			if !bias2.is_static
				update_bias!(bias2, field2.CV)
			end
		end

		return accept_swap
	end 

	function swap!(a::Gaugefield, b::Gaugefield)
		NX, NT, _ = size(a)

		for μ in 1:2
			for it in 1:NT
				for ix in 1:NX
					a[ix,it,μ], b[ix,it,μ] = b[ix,it,μ], a[ix,it,μ]
				end
			end
		end

		return nothing
	end

end