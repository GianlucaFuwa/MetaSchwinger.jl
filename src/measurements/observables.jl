module Observables
    using Arpack

    import ..Gaugefields: circshift, Gaugefield, plaquette
    import ..Gaugefields: wilson_loop

    function meta_charge(g::Gaugefield)
		NX, NT, _ = size(g)
		q = 0.0

		for it in 1:NT
			for ix in 1:NX
			    q += sin(plaquette(g, ix, it))
			end
		end

		return q / 2π
	end

    function topological_charge(g::Gaugefield)
        NX, NT, _ = size(g)
        q = 0.0

        for it in 1:NT
            for ix in 1:NX
                @inbounds q += imag(log(cis(plaquette(g, ix, it))))
            end
        end

        return round(Int, q / 2π, RoundNearestTiesAway)
    end

    function wilson_loop_all(g::Gaugefield, LT)
        NX, _, _ = size(g)
        wils = zeros(ComplexF64, NX)

        for LX in 1:NX
            wils[LX] += wilson_loop_avg(g, LX, LT, tempgauge = true)
        end

        return real.(wils), imag.(wils)
    end

    function poly_loop_avg(g::Gaugefield)
        NX, _, _ = size(g)
        poly_re = 0.0
        poly_im = 0.0

        for ix in 1:NX
            pre, pim = poly_loop(g, ix)
            poly_re += pre
            poly_im += pim
        end

        return poly_re / NX, poly_im / NX
    end

    function poly_loop(g::Gaugefield, ix)
		_, NT, _ = size(g)
		poly = 0.0

		for it in 1:NT
			poly += g[ix,it,2]
		end
		
		return cos(poly), sin(poly)
	end 

	function wilson_loop_avg(g::Gaugefield, LX, LT; tempgauge = false)
		NX, NT, _ = size(g)
		wils = zeros(NX, NT)
		x_up_side = zeros(NX, NT)
		t_up_side = zeros(NX, NT)
		x_down_side = zeros(NX, NT)
		t_down_side = zeros(NX, NT)

		if tempgauge == true

			for i in 0:LX-1
				x_up_side += circshift(g, (i,0), 1)
				x_down_side -= circshift(g, (i,LT), 1)
			end

			wils = x_up_side + x_down_side
			return sum(cis.(wils)) / g.NV
		elseif tempgauge == false

			for i in 0:LX-1
				x_up_side += circshift(g, (i,0), 1)
				x_down_side -= circshift(g, (i,LT), 1)
			end
			
			for i in 0:LT-1
				t_up_side += circshift(g, (LX,i), 1)
				t_down_side -= circshift(g, (0,i), 1)
			end

			wils = x_up_side + t_up_side + x_down_side + t_down_side
			return sum(cis.(wils)) / g.NV
		else
			error("tempgauge in Wilson loop calculation can only either be true or false")
		end
	end

end