function tempering_heatbath!(fields, biases, rng)
    numinstances = length(fields)

    accepted_swaps = 0
    for i in 2:numinstances
        accepted_swaps += tempering_swap!(fields[1], fields[i], biases[1], biases[i], rng)
    end

    return accepted_swaps
end

function tempering_swap!(field1, field2, bias1, bias2, rng)
	CV1 = field1.CV
	CV2 = field2.CV

	ΔCV = CV1 - CV2
	ΔSg = field1.Sg - field2.Sg
	ΔV1 = bias1(CV2) - bias1(CV1)
	ΔV2 = bias2(CV1) - bias2(CV2)
	accept_swap = rand(rng) ≤ exp(-ΔV1-ΔV2)

	if accept_swap
		swap!(field1, field2)
		field2.CV += ΔCV
		field1.CV -= ΔCV
		field2.Sg += ΔSg
		field1.Sg -= ΔSg

        update_bias!(bias1, field1.CV)
        update_bias!(bias2, field2.CV)
	end

	return accept_swap
end

function swap!(a, b)
	NX, NT = size(a)

    for it in 1:NT
        for ix in 1:NX
            for μ in 1:2
                a[μ,ix,it], b[μ,ix,it] = b[μ,ix,it], a[μ,ix,it]
			end
		end
	end

	return nothing
end
