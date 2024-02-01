module Observables
    using Arpack

    import ..Gaugefields: Gaugefield, plaquette

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

    function wilson_loop(u::Gaugefield, LT, LX)
        NX, NT = size(u)
        temp = zeros(NX, NT)
        wl = wilson_loop!(temp, u, LT, LX)
        return wl
    end

    function wilson_loop!(temp, u::Gaugefield, LX, LT)
        for ix in 0:LX-1
            temp .+= circshift(view(u.U, 1, :, :), (-ix, 0))
            temp .-= circshift(view(u.U, 1, :, :), (-ix, -LT))
        end

        for it in 0:LT-1
            temp .+= circshift(view(u.U, 2, :, :), (-LX, -it))
            temp .-= circshift(view(u.U, 2, :, :), (0, -it))
        end

        return sum(cos.(temp)) / u.NV
    end
end
