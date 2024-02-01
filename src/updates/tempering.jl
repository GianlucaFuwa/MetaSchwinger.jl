function tempering_heatbath!(Us, biases, rng)
    numinstances = length(Us)

    accepted_swaps = 0
    for i in 2:numinstances
        accepted_swaps += tempering_swap!(Us[1], Us[i], biases[1], biases[i], rng)
    end

    return accepted_swaps
end

function tempering_swap!(U1, U2, bias1, bias2, rng)
	CV1 = U1.CV
	CV2 = U2.CV

	ΔCV = CV1 - CV2
	ΔSg = U1.Sg - U2.Sg
	ΔV1 = bias1(CV2) - bias1(CV1)
	ΔV2 = bias2(CV1) - bias2(CV2)
	accept_swap = rand(rng) ≤ exp(-ΔV1-ΔV2)

	if accept_swap
		swap!(U1, U2)
		U2.CV += ΔCV
		U1.CV -= ΔCV
		U2.Sg += ΔSg
		U1.Sg -= ΔSg

        update_bias!(bias1, U1.CV)
        update_bias!(bias2, U2.CV)
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
