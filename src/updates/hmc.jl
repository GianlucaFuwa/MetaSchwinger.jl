abstract type AbstractIntegrator end

struct HMCUpdate{TI} <: AbstractUpdate
    steps::Int64
    Δτ::Float64
    ϕ::Float64
    P::Liefield
    _temp_U::Gaugefield
    _temp_force::Liefield

    function HMCUpdate(U, integrator, steps, Δτ, ϕ)
        P = Liefield(U)
        _temp_U = similar(U)
        _temp_force = Liefield(U)

        TI = getfield(Updates, Symbol(integrator))
        return new{TI}(steps, Δτ, ϕ, P, _temp_U, _temp_force)
    end
end

include("hmc_integrators.jl")

function update!(updatemethod::HMCUpdate{TI}, U, rng; bias=nothing, metro_test=true) where {TI}
    U_old = updatemethod._temp_U
    substitute_U!(U_old, U)
    gaussian_momenta!(updatemethod.P, updatemethod.ϕ, rng)

    Sg_old = U.Sg
    CV_old = U.CV
    trP²_old = calc_kinetic_energy(updatemethod.P)

    evolve!(TI(), U, updatemethod, bias)

    Sg_new = calc_Sg(U)
    trP²_new = calc_kinetic_energy(updatemethod.P)

    ΔSg = Sg_new - Sg_old
    ΔP² = trP²_new - trP²_old

    if bias === nothing
        CV_new = CV_old
        ΔV = 0
    else
        CV_new = calc_CV(U)
        ΔV = 19.5cos(1.044π*CV_new) - 19.5cos(1.044π*CV_old)#bias(CV_new) - bias(CV_old)
    end

    ΔH = ΔP² + ΔSg + ΔV
    println("ΔP² = ", ΔP²)
    println("ΔSg = ", ΔSg)
    println("ΔV = ", ΔV)
    println("ΔH = ", ΔH)

    accept = metro_test ? rand(rng)<=exp(-ΔH) : true
    if accept
        println("accepted")
        U.Sg = Sg_new
        U.CV = CV_new
    else
        println("rejected")
        substitute_U!(U, U_old)
        updatemethod.ϕ!=π/2 && flip_momenta!(updatemethod.P)
    end

    return accept
end

function updateU!(U, method, fac)
    ϵ = method.Δτ * fac
    P = method.P
    add!(U, P, ϵ)
    return nothing
end

function updateP!(U, method, fac, bias)
    ϵ = method.Δτ * fac
    P = method.P
    force = method._temp_force

    if bias !== nothing
        calc_dVdU!(force, U, bias)
        add!(P, force, ϵ)
    end

    calc_dSdU!(force, U)
    add!(P, force, ϵ)
    return nothing
end

function calc_dSdU!(dSdU, U)
    NX, NT = size(U)
    prefactor = -U.β

    for it in 1:NT
        for ix in 1:NX
            for μ in 1:2
                A = staple(U, μ, ix, it)
                UA′ = cis(U[μ,ix,it]) * A'
                dSdU[μ,ix,it] = prefactor * imag(UA′)
            end
        end
    end

    return nothing
end

function calc_dVdU!(dVdU, U, bias)
    NX, NT = size(U)
    prefactor = 1/2π
    Q = calc_CV(U)
    bias_der = -19.5*1.044π * sin(1.044π*Q)#∂V∂Q(bias, Q)

    for it in 1:NT
        for ix in 1:NX
            for μ in 1:2
                A = staple(U, μ, ix, it)
                UA′ = cis(U[μ,ix,it]) * A'
                dVdU[μ,ix,it] = bias_der * prefactor * sin(U[μ,ix,it]) * 0
            end
        end
    end

    return nothing
end
