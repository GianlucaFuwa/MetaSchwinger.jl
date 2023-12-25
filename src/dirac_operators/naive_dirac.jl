struct NaiveDiracOperator <: AbstractDiracOperator
    Dop::Matrix{ComplexF64}
    expiU::Array{ComplexF64, 3}
    mass::Float64
    BC::Vector{Int64}

    function NaiveDiracOperator(N, mass, BC)
        NX, NT = N
        len = NX * NT * 2
        Dop = zeros(ComplexF64, len, len)
        expiU = Array{ComplexF64, 3}(undef, NX, NT, 2)
        return new(Dop, expiU, mass, BC)
    end
end

function (D::NaiveDiracOperator)(U::Gaugefield)
    NX, NT = size(U)
    mass = D.mass
    BC_x, BC_t = D.BC
    expiU = D.expiU
    Dop = D.Dop

    clear!(D)
    expiU .= cis(U)

    @batch for nt in 1:NT
        for nx in 1:NX
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

            view(Dop, n:(n+1), n:(n+1)) .+= mass * eye2
            view(Dop, n:(n+1), npx:(npx+1)) .-= 0.5 * γ1 * expiU[nx,nt,1] * BC_xp
            view(Dop, n:(n+1), nmx:(nmx+1)) .-= 0.5 * γ1 * conj(expiU[nx_min,nt,1]) * BC_xm
            view(Dop, n:(n+1), npt:(npt+1)) .-= 0.5 * γ2 * expiU[nx,nt,2] * BC_tp
            view(Dop, n:(n+1), nmt:(nmt+1)) .-= 0.5 * γ2 * conj(expiU[nx,nt_min,2]) * BC_tm
        end
    end

    return nothing
end
