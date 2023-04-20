module Observables

    import ..Gaugefields: Gaugefield,plaquette,staple,wilson_loop,wilson_loop_avg,poly_loop

    function MetaCharge(g::Gaugefield)
		NX,NT,_ = size(g)
		q = 0.
		for it=1:NT
			for ix=1:NX
			q += sin(plaquette(g,ix,it))
			end
		end
		return q / 2pi
	end

    function TopCharge(g::Gaugefield)
        NX,NT,_ = size(g)
        q = 0.0 + 0.0im
        for it=1:NT
            for ix=1:NX
                @inbounds q += log(exp(plaquette(g,ix,it)im))
            end
        end
        return round(Int,imag(q)/2pi,RoundNearestTiesAway)
    end

    function wilson_loop_all(g::Gaugefield,LT::Int64)
        NX,_,_ = size(g)
        wils = zeros(ComplexF64,NX)
        for LX=1:NX
            wils[LX] += wilson_loop_avg(g,LX,LT,tempgauge=true)
        end
        return real.(wils),imag.(wils)
    end

    function poly_loop_avg(g::Gaugefield)
        NX,_,_ = size(g)
        poly_re = 0.0
        poly_im = 0.0
        for ix=1:NX
            pre,pim = poly_loop(g,ix)
            poly_re += pre
            poly_im += pim
        end 
        return poly_re/NX,poly_im/NX
    end

end