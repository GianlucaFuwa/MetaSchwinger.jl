function metropolis!(field::Gaugefield, ix, it, μ, rng, ϵ, multi_hit)
	accept = 0

	for _ in 1:multi_hit
		dU = randn(rng) * ϵ
		ΔSg = local_action_diff(field, ix, it, μ, dU)
		ΔCV = local_metacharge_diff(field, ix, it, μ, dU)

		if rand(rng) ≤ exp(-ΔSg)
			field[ix,it,μ] += dU
			field.Sg += ΔSg
			field.CV += ΔCV
			accept += 1
		end
	end

	return accept
end

function metropolis_meta!(field::Gaugefield, bias::BiasPotential, ix, it, μ, rng, ϵ, multi_hit)
	accept = 0

	for _ in 1:multi_hit
		dU = randn(rng) * ϵ
		ΔSg = local_action_diff(field, ix, it, μ, dU)
		ΔCV = local_metacharge_diff(field, ix, it, μ, dU)

		oldCV = field.CV
		newCV = field.CV + ΔCV

		ΔV = bias(newCV) - bias(oldCV)

		if rand(rng) ≤ exp(-ΔSg - ΔV)
			field[ix,it,μ] += dU
			field.Sg += ΔSg
			field.CV += ΔCV
			accept += 1
		end
	end

	return accept
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

function sweep_meta!(field::Gaugefield, bias::BiasPotential, rng, ϵ, multi_hit)
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

function adjusted_ϵ(ϵ, numaccepts, metro_norm, metro_target_acc)
	return 	ϵ + (numaccepts * metro_norm - metro_target_acc) * 0.2
end

function local_action_diff(g::Gaugefield, ix, it, μ, dU)
	NX, NT = size(g)
	it_min = mod1(it - 1, NT)
	it_plu = mod1(it + 1, NT)
	ix_min = mod1(ix - 1, NX)
	ix_plu = mod1(ix + 1, NX)

	if μ == 1
		a = g[ix,it    ,1]+dU + g[ix_plu,it    ,2] - g[ix,it_plu,1]    - g[ix,it    ,2]
		b = g[ix,it_min,1]    + g[ix_plu,it_min,2] - g[ix,it    ,1]-dU - g[ix,it_min,2]
		c = g[ix,it    ,1]    + g[ix_plu,it    ,2] - g[ix,it_plu,1]    - g[ix,it    ,2]
		d = g[ix,it_min,1]    + g[ix_plu,it_min,2] - g[ix,it    ,1]    - g[ix,it_min,2]
	elseif μ == 2
		a = g[ix_min,it,1] + g[ix    ,it,2]+dU - g[ix_min,it_plu,1] - g[ix_min,it,2]
		b = g[ix    ,it,1] + g[ix_plu,it,2]    - g[ix    ,it_plu,1] - g[ix    ,it,2]-dU
		c = g[ix_min,it,1] + g[ix    ,it,2]    - g[ix_min,it_plu,1] - g[ix_min,it,2]
		d = g[ix    ,it,1] + g[ix_plu,it,2]    - g[ix    ,it_plu,1] - g[ix    ,it,2]
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
		a = g[ix,it    ,1]+dU + g[ix_plu,it    ,2] - g[ix,it_plu,1]    - g[ix,it    ,2]
		b = g[ix,it_min,1]    + g[ix_plu,it_min,2] - g[ix,it    ,1]-dU - g[ix,it_min,2]
		c = g[ix,it    ,1]    + g[ix_plu,it    ,2] - g[ix,it_plu,1]    - g[ix,it    ,2]
		d = g[ix,it_min,1]    + g[ix_plu,it_min,2] - g[ix,it    ,1]    - g[ix,it_min,2]
	elseif μ == 2
		a = g[ix_min,it,1] + g[ix    ,it,2]+dU - g[ix_min,it_plu,1] - g[ix_min,it,2]
		b = g[ix    ,it,1] + g[ix_plu,it,2]    - g[ix    ,it_plu,1] - g[ix    ,it,2]-dU
		c = g[ix_min,it,1] + g[ix    ,it,2]    - g[ix_min,it_plu,1] - g[ix_min,it,2]
		d = g[ix    ,it,1] + g[ix_plu,it,2]    - g[ix    ,it_plu,1] - g[ix    ,it,2]
	end

	return 1 / 2π * (sin(a) + sin(b) - sin(c) - sin(d))
end
