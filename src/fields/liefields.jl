struct Liefield  <: Abstractfield
    NX::Int64
    NT::Int64
    NV::Int64
    U::Array{Float64,3}

    function Liefield(NX, NT)
        NV = NX * NT
        U = zeros(Float64, 2, NX, NT)
        return new(NX, NT, NV, U)
    end
end

Liefield(U) = Liefield(size(U)...)

function gaussian_momenta!(p::Liefield, ϕ, rng)
    NX, NT = size(p)
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)

    for it in 1:NT
        for ix in 1:NX
            for μ in 1:2
                p[μ,ix,it] = cosϕ*p[μ,ix,it] + sinϕ*randn(rng)
            end
        end
    end

    return nothing
end

function calc_kinetic_energy(p::Liefield)
    NX, NT = size(p)
    P = p.U

    kin = 0.0
    @turbo for it in 1:NT
        for ix in 1:NX
            for μ in 1:2
                kin += P[μ,ix,it]^2
            end
        end
    end

    return 0.5kin
end
