module Observables
    using Arpack

    import ..Gaugefields: Gaugefield, plaquette, wilson_loop

    function topological_charge(u::Gaugefield)
        NX, NT = size(u)
        q = 0.0

        for it in 1:NT
            for ix in 1:NX
                q += imag(log(cis(plaquette(u, ix, it))))
            end
        end

        return round(Int, q / 2π, RoundNearestTiesAway)
    end

    function wilson_loop_all(u::Gaugefield, LT)
        NX, _ = size(u)
        wils = zeros(ComplexF64, NX)

        for LX in 1:NX
            wils[LX] += wilson_loop_avg(u, LX, LT, tempgauge=true)
        end

        return real.(wils), imag.(wils)
    end

    function poly_loop_avg(u::Gaugefield)
        NX, _ = size(u)
        poly_re = 0.0
        poly_im = 0.0

        for ix in 1:NX
            pre, pim = poly_loop(u, ix)
            poly_re += pre
            poly_im += pim
        end

        return poly_re/NX, poly_im/NX
    end

    function poly_loop(u::Gaugefield, ix)
		_, NT = size(u)
		poly = 0.0

		for it in 1:NT
			poly += u[2,ix,it]
		end

		return cos(poly), sin(poly)
	end

	function wilson_loop_avg(u::Gaugefield, LX, LT; tempgauge=false)
		NX, NT = size(u)
		wils = zeros(NX, NT)
		x_up_side = zeros(NX, NT)
		t_up_side = zeros(NX, NT)
		x_down_side = zeros(NX, NT)
		t_down_side = zeros(NX, NT)

		if tempgauge == true

			for i in 0:LX-1
				x_up_side += circshift(u, (i,0), 1)
				x_down_side -= circshift(u, (i,LT), 1)
			end

			wils = x_up_side + x_down_side
			return sum(cis.(wils)) / u.NV
		elseif tempgauge == false

			for i in 0:LX-1
				x_up_side += circshift(u, (i,0), 1)
				x_down_side -= circshift(u, (i,LT), 1)
			end

			for i in 0:LT-1
				t_up_side += circshift(u, (LX,i), 1)
				t_down_side -= circshift(u, (0,i), 1)
			end

			wils = x_up_side + t_up_side + x_down_side + t_down_side
			return sum(cis.(wils)) / u.NV
		else
			error("tempgauge in Wilson loop calculation can only either be true or false")
		end
	end

end
