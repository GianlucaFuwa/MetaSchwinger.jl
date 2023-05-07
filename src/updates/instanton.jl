module InstantonUpdate
    using Random

    import ..Gaugefields: calc_Sg, Gaugefield, substitute_U!
	#import ..Verbose_print: println_verbose, Verbose_

    function instanton_update!(U::Gaugefield, Q, rng)
        NX, NT, _ = size(U)
        Uold = deepcopy(U.U)

        Sg_old = U.Sg
    
        for it in 1:NX
            for ix in 1:NT
                U[ix,it,1] -= Q * it * 2π / (NX * NT)
                U[ix,it,2] += Q * ix * (it == NT) * 2π / NX
            end
        end

        Sg_new = calc_Sg(U)
        ΔSg = Sg_new - Sg_old
        accept = rand(rng) ≤ exp(-ΔSg)
        
        if accept
            U.Sg += ΔSg
        else
            substitute_U!(U, Uold)
        end

        return accept
    end

end
        