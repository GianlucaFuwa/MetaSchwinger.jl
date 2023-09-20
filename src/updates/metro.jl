mutable struct MetroUpdate <: AbstractUpdate
    ϵ::Float64
    multi_hit::Int64
    target_acc::Float64
    norm::Float64

    function MetroUpdate(U, ϵ, multi_hit, target_acc)
        norm = 1 / (U.NV * 2 * multi_hit)
        return new(ϵ, multi_hit, target_acc, norm)
    end
end

function update!(updatemethod::MetroUpdate, U, rng; bias=nothing, kwargs...)
    ϵ = updatemethod.ϵ
    multi_hit = updatemethod.multi_hit
    target_acc = updatemethod.target_acc
    norm = updatemethod.norm

    if bias === nothing
        numaccepts = sweep!(U, rng, ϵ, multi_hit)
    else
        numaccepts = sweep_meta!(U, bias, rng, ϵ, multi_hit)
    end

    updatemethod.ϵ = adjusted_ϵ(ϵ, numaccepts, norm, target_acc)
    return numaccepts*norm
end

function sweep!(field::Gaugefield, rng, ϵ, multi_hit)
	NX, NT = size(field)
	numaccepts = 0

	for eo in 0:1
		for ix in (eo+1):2:NX
			for it in 1:NT
				numaccepts += metropolis!(field, ix, it, 2, rng, ϵ, multi_hit)
			end
		end
	end

	for eo in 0:1
		for it in (eo+1):2:NT
			for ix in 1:NX
				numaccepts += metropolis!(field, ix, it, 1, rng, ϵ, multi_hit)
			end
		end
	end

	return numaccepts
end

function sweep_meta!(field::Gaugefield, bias, rng, ϵ, multi_hit)
	NX, NT = size(field)
	numaccepts = 0

	for eo in 0:1
		for ix in (eo+1):2:NX
			for it in 1:NT
				numaccepts += metropolis_meta!(field, bias, ix, it, 2, rng, ϵ, multi_hit)
			end
		end
	end

	for eo in 0:1
		for it in (eo+1):2:NT
			for ix in 1:NX
				numaccepts += metropolis_meta!(field, bias, ix, it, 1, rng, ϵ, multi_hit)
			end
		end
	end

	return numaccepts
end

function metropolis!(field::Gaugefield, ix, it, μ, rng, ϵ, multi_hit)
	accept = 0

	for _ in 1:multi_hit
		dU = randn(rng) * ϵ
		ΔSg = local_action_diff(field, ix, it, μ, dU)
		ΔCV = local_metacharge_diff(field, ix, it, μ, dU)

		if rand(rng) ≤ exp(-ΔSg)
			field[μ,ix,it] += dU
			field.Sg += ΔSg
			field.CV += ΔCV
			accept += 1
		end
	end

	return accept
end

function metropolis_meta!(field::Gaugefield, bias, ix, it, μ, rng, ϵ, multi_hit)
	accept = 0

	for _ in 1:multi_hit
		dU = randn(rng) * ϵ
		ΔSg = local_action_diff(field, ix, it, μ, dU)
		ΔCV = local_metacharge_diff(field, ix, it, μ, dU)

		oldCV = field.CV
		newCV = field.CV + ΔCV

		ΔV = bias(newCV) - bias(oldCV)

		if rand(rng) ≤ exp(-ΔSg - ΔV)
			field[μ,ix,it] += dU
			field.Sg += ΔSg
			field.CV += ΔCV
			accept += 1
		end
	end

	return accept
end

function adjusted_ϵ(ϵ, numaccepts, metro_norm, metro_target_acc)
	return 	mod(ϵ + (numaccepts*metro_norm - metro_target_acc) * 0.2, 2π)
end

function local_action_diff(g::Gaugefield, ix, it, μ, dU)
	NX, NT = size(g)
	it_min = mod1(it - 1, NT)
	it_plu = mod1(it + 1, NT)
	ix_min = mod1(ix - 1, NX)
	ix_plu = mod1(ix + 1, NX)

	if μ == 1
		a = g[1,ix,it    ]+dU + g[2,ix_plu,it    ] - g[1,ix,it_plu]    - g[2,ix,it    ]
		b = g[1,ix,it_min]    + g[2,ix_plu,it_min] - g[1,ix,it    ]-dU - g[2,ix,it_min]
		c = g[1,ix,it    ]    + g[2,ix_plu,it    ] - g[1,ix,it_plu]    - g[2,ix,it    ]
		d = g[1,ix,it_min]    + g[2,ix_plu,it_min] - g[1,ix,it    ]    - g[2,ix,it_min]
	elseif μ == 2
		a = g[1,ix_min,it] + g[2,ix    ,it]+dU - g[1,ix_min,it_plu] - g[2,ix_min,it]
		b = g[1,ix    ,it] + g[2,ix_plu,it]    - g[1,ix    ,it_plu] - g[2,ix    ,it]-dU
		c = g[1,ix_min,it] + g[2,ix    ,it]    - g[1,ix_min,it_plu] - g[2,ix_min,it]
		d = g[1,ix    ,it] + g[2,ix_plu,it]    - g[1,ix    ,it_plu] - g[2,ix    ,it]
	end

	return -g.β * (cos(a) + cos(b) - cos(c) - cos(d))
end

function local_metacharge_diff(g::Gaugefield, ix, it, μ, dU)
	NX, NT = size(g)
	it_min = mod1(it - 1, NT)
	it_plu = mod1(it + 1, NT)
	ix_min = mod1(ix - 1, NX)
	ix_plu = mod1(ix + 1, NX)

	if μ == 1
		a = g[1,ix,it    ]+dU + g[2,ix_plu,it    ] - g[1,ix,it_plu]    - g[2,ix,it    ]
		b = g[1,ix,it_min]    + g[2,ix_plu,it_min] - g[1,ix,it    ]-dU - g[2,ix,it_min]
		c = g[1,ix,it    ]    + g[2,ix_plu,it    ] - g[1,ix,it_plu]    - g[2,ix,it    ]
		d = g[1,ix,it_min]    + g[2,ix_plu,it_min] - g[1,ix,it    ]    - g[2,ix,it_min]
	elseif μ == 2
		a = g[1,ix_min,it] + g[2,ix    ,it]+dU - g[1,ix_min,it_plu] - g[2,ix_min,it]
		b = g[1,ix    ,it] + g[2,ix_plu,it]    - g[1,ix    ,it_plu] - g[2,ix    ,it]-dU
		c = g[1,ix_min,it] + g[2,ix    ,it]    - g[1,ix_min,it_plu] - g[2,ix_min,it]
		d = g[1,ix    ,it] + g[2,ix_plu,it]    - g[1,ix    ,it_plu] - g[2,ix    ,it]
	end

	return 1 / 2π * (sin(a) + sin(b) - sin(c) - sin(d))
end
