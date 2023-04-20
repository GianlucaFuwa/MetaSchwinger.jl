module Tempering
    using Random
    
    import ..Gaugefields: Gaugefield,swap!,add_Sg!,add_CV!
    import ..Metadynamics: Bias_potential,update_bias!,DeltaV

    function tempering_swap!(field0::Gaugefield, field1::Gaugefield, bias::Bias_potential, rng::Xoshiro)
		CV0 = field0.CV
		CV1 = field1.CV

		ΔCV = CV0 - CV1
		ΔSg = field0.Sg - field1.Sg
		ΔV = DeltaV(bias, CV0, CV1)
		accept_swap = rand(rng) ≤ exp(ΔV)
		if accept_swap
			swap!(field0, field1)
			add_CV!(field1, ΔCV)
			add_CV!(field0, -ΔCV)
			add_Sg!(field1, ΔSg)
			add_Sg!(field0, -ΔSg)
			if !bias.is_static
				update_bias!(bias, field1.CV)
			end
		end
		return accept_swap
	end 

end