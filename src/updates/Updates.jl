module Updates
    using LoopVectorization
    using Random

    import ..BiasModule: Bias, ∂V∂Q
    import ..Gaugefields: Gaugefield, Liefield, add!, calc_kinetic_energy, calc_CV, calc_Sg,
        gaussian_momenta!, recalc_CV!, substitute_U!, staple
    import ..Parameters: ParameterSet

    abstract type AbstractUpdate end

    update!(::AbstractUpdate, args...) = nothing
    update!(::Nothing, args...) = nothing

    include("metro.jl")
    include("hmc.jl")
    include("instanton.jl")
    include("tempering.jl")

    function Updatemethod(p::ParameterSet, U)
        updatemethod = Updatemethod(
            U,
            p.update_method,
            p.metro_ϵ,
            p.metro_multi_hit,
            p.metro_target_acc,
            p.hmc_integrator,
            p.hmc_steps,
            p.hmc_Δτ,
            p.hmc_ϕ,
        )
        return updatemethod
    end

    function Updatemethod(
        U,
        update_method,
        metro_ϵ = 0.2,
        metro_multi_hit = 1,
        metro_target_acc = 0.5,
        hmc_integrator = "OMF4",
        hmc_steps = 10,
        hmc_Δτ = 0.1,
        hmc_ϕ = π/2,
    )
        if update_method ∈ ["hmc", "HMC"]
            updatemethod = HMCUpdate(U, hmc_integrator, hmc_steps, hmc_Δτ, hmc_ϕ)
        elseif update_method ∈ ["metro", "Metro", "metropolis", "Metropolis"]
            updatemethod = MetroUpdate(U, metro_ϵ, metro_multi_hit, metro_target_acc)
        else
            error("update method $(update_method) not supported")
        end

        return updatemethod
    end
end
