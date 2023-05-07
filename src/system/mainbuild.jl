module Mainbuild
    using Dates
    using DelimitedFiles
    using InteractiveUtils
    
    import ..DiracOperatorModule: dirac_operator
    import ..Gaugefields: Gaugefield, recalc_Sg!, recalc_CV!
    import ..MetroUpdate: adjusted_ϵ, sweep!, sweep_meta!
    import ..InstantonUpdate: instanton_update!
    import ..Metadynamics: BiasPotential, parametric_to_bias!, update_bias!, write_to_file!
    import ..Measurements: calc_weights, MeasurementSet, measurements
    import ..System_parameters: parameterloading, Params, ParamSet
    import ..System_parameters: physical, meta, param_meta, sim, mc, dirac, meas, system
    import ..Tempering: tempering_swap!
    import ..Verbose_print: print2file, println_verbose, Verbose_

    function run_build(filenamein::String)
        filename = filenamein
        include(abspath(filename))
        params_set = ParamSet(physical, meta, param_meta, sim, mc, dirac, meas, system)

        run_build(params_set)

        return nothing
    end

    function run_build(params_set::ParamSet)
        params = parameterloading(params_set)

        field = Gaugefield(params)
        bias = BiasPotential(params)
        run_build!(field,bias,params)

        return nothing
    end

    function run_build!(field::Gaugefield, bias::BiasPotential, params::Params)
        verbose = Verbose_(params.logfile)
        println_verbose(verbose, "# ", pwd())
        println_verbose(verbose, "# ", Dates.now())
        versioninfo(verbose)

        measset = MeasurementSet(params.measure_dir, meas_calls = params.meas_calls)
        ϵ = params.ϵ_metro
        multi_hit = params.multi_hit
        metro_target_acc = params.metro_target_acc
        metro_norm = 1 / (field.NV * 2 * multi_hit)
        rng = params.randomseeds[1]
        Nsweeps = params.Nsweeps

        if params.initial == "hot"
            field.U = rand(size(field.U) .- 0.5) * 2 * 2π
            recalc_Sg!(field)
            recalc_CV!(field)
        end 

        for therm in 1:params.Ntherm
            tmp = sweep!(field, rng, ϵ, multi_hit)
            ϵ = adjusted_ϵ(ϵ, tmp, metro_norm, metro_target_acc)
        end
        
        if bias.parametric == true
            parametric_to_bias!(bias)
        end
        
        numaccepts = 0

        for itrj in 1:Nsweeps
            tmp = sweep_meta!(field, bias, rng, ϵ, multi_hit)
            ϵ = adjusted_ϵ(ϵ, tmp, metro_norm, metro_target_acc)
            numaccepts += tmp
            update_bias!(bias, field.CV, itrj = itrj)

            if params.veryverbose
                println_verbose(
                    verbose,
                    " ",
                    itrj,
                    " ",
                    100 * numaccepts / trj / field.NV / 2,
                    "%",
                    " # trj accrate"
                )
            end

            measurements(itrj, field, measset)
        end

        println_verbose(verbose, "Final ϵ_metro: ", ϵ)
        println_verbose(
            verbose,
            "Acceptance rate: ",
            100 * numaccepts * metro_norm / Nsweeps,
            "%"
        )
        println_verbose(verbose, "Exceeded Count: ", bias.exceeded_count)

        write_to_file!(bias)
        flush(bias)

        flush(stdout)
        flush(verbose)
        return nothing
    end

end
