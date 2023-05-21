struct WilsonDiracOperator <: AbstractDiracOperator
    mass::Float64
    BC::Vector{Int64}
end

function (D::WilsonDiracOperator)(U::Gaugefield)
    NX, NT, d = size(U)
    mass = D.mass
    BC_x, BC_t = D.BC

    len = NX * NT * d
    expiU = cis(U)

    Dwils = zeros(ComplexF64, len, len)

    for nx in 1:NX
        for nt in 1:NT
            n = flat_index(NX, nx, nt)

            BC_xp = ifelse(nx == NX, BC_x, 1)
            BC_xm = ifelse(nx == 1, BC_x, 1)
            BC_tp = ifelse(nt == NT, BC_t, 1)
            BC_tm = ifelse(nt == 1, BC_t, 1)

            nx_min = mod1(nx - 1, NX)
            nt_min = mod1(nt - 1, NT)

            npx = flat_index(NX, mod1(nx + 1, NX), nt)
            nmx = flat_index(NX, nx_min, nt)
            npt = flat_index(NX, nx, mod1(nt + 1, NT))
            nmt = flat_index(NX, nx, nt_min)
            
            view(Dwils, n:(n+1), n:(n+1)) .+= (mass + 2) * eye2
            view(Dwils, n:(n+1), npx:(npx+1)) .-= 0.5 * (eye2 - γ1) * expiU[nx,nt,1] * BC_xp
            view(Dwils, n:(n+1), nmx:(nmx+1)) .-= 0.5 * (eye2 + γ1) * conj(expiU[nx_min,nt,1]) * BC_xm
            view(Dwils, n:(n+1), npt:(npt+1)) .-= 0.5 * (eye2 - γ2) * expiU[nx,nt,2] * BC_tp
            view(Dwils, n:(n+1), nmt:(nmt+1)) .-= 0.5 * (eye2 + γ2) * conj(expiU[nx,nt_min,2]) * BC_tm
        end
    end

    return Dwils
end


