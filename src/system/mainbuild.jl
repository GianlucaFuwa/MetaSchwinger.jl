module Mainbuild
    using Base.Threads
    using Dates
    using DelimitedFiles
    using InteractiveUtils
    
    import ..DiracOperators: dirac_operator
    import ..Gaugefields: Gaugefield, recalc_Sg!, recalc_CV!
    import ..Updates: adjusted_ϵ, sweep!, sweep_meta!
    import ..Updates: instanton!, instanton_update!, tempering_swap!
    import ..Metadynamics: BiasPotential, clear!, parametric_to_bias!, update_bias!, write_to_file!
    import ..Measurements: calc_weights, MeasurementSet, measurements
    import ..Parameters: parameterloading, ParameterSet, ParamSet
    import ..Parameters: physical, meta, param_meta, update, dirac, meas, system
    import ..VerbosePrint: print2file, println_verbose, Verbose_

    function run_build(filenamein::String)
        filename = filenamein
        include(abspath(filename))
        params_set = ParamSet(physical, meta, param_meta, update, dirac, meas, system)
        run_build(params_set)
        return nothing
    end

    function run_build(params_set::ParamSet)
        params = parameterloading(params_set)

        if params.tempering_enabled
            fields = Vector{Gaugefield}(undef, 0)
            biases = Vector{BiasPotential}(undef, 0)

            for _ in 1:params.numinstances
                push!(fields, Gaugefield(params))
            end

            for _ in 1:params.numinstances
                push!(
                    biases,
                    BiasPotential(params, instance = 1),
                )
            end
            
            cum_bias = BiasPotential(params)
            run_tempered_build!(fields, biases, cum_bias, params)
        elseif !params.tempering_enabled
            field = Gaugefield(params)
            bias = BiasPotential(params)
            run_build!(field, bias, params)
        end

        return nothing
    end

    function run_build!(field::Gaugefield, bias::BiasPotential, params::ParameterSet)
        verbose = Verbose_(params.logfile)
        println_verbose(verbose, "# ", pwd())
        println_verbose(verbose, "# ", Dates.now())
        versioninfo(verbose)

        measset = MeasurementSet(params.measure_dir, meas_calls = params.meas_calls)
        ϵ = params.ϵ_metro
        multi_hit = params.metro_multi_hit
        metro_target_acc = params.metro_target_acc
        metro_norm = 1 / (field.NV * 2 * multi_hit)
        rng = params.randomseeds[1]
        Nsweeps = params.Nsweeps

        if params.initial == "hot"
            field.U = rand(size(field.U) .- 0.5) * 2 * 2π
            recalc_Sg!(field)
            recalc_CV!(field)
        end 

        for _ in 1:params.Ntherm
            tmp = sweep!(field, rng, ϵ, multi_hit)
            ϵ = adjusted_ϵ(ϵ, tmp, metro_norm, metro_target_acc)
        end
        
        if bias.parametric == true
            #parametric_to_bias!(bias)
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
                    100 * numaccepts / itrj / field.NV / 2,
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

        write_to_file!(bias, bias.fp)
        flush(bias)

        flush(stdout)
        flush(verbose)
        return nothing
    end

    function run_tempered_build!(
        fields::Vector{Gaugefield},
        biases::Vector{BiasPotential},
        cum_bias::BiasPotential,
        params::ParameterSet,
    )
        numinstances = params.numinstances
        verbose = Verbose_(params.logfile)
        println_verbose(verbose, "# ", pwd())
        println_verbose(verbose, "# ", Dates.now())
        versioninfo(verbose)
        println_verbose(
            verbose,
            "Parallel tempered run with ",
            numinstances,
            " MetaD instances",
        )

        meassets = Vector{MeasurementSet}(undef, 0)

        for i in 1:numinstances
            push!(
                meassets,
                MeasurementSet(params.measure_dir, meas_calls = params.meas_calls,
                instance = "_$(i)")
            )
        end

        ϵ = params.ϵ_metro
        multi_hit = params.metro_multi_hit
        metro_target_acc = params.metro_target_acc
        metro_norm = 1 / (fields[1].NV * 2 * multi_hit)
        rng = params.randomseeds

        if params.initial == "hot"
            for i in 1:numinstances
            fields[i].U = rand(size(fields[i].U) .- 0.5) * 2 * 2π
            recalc_Sg!(fields[i])
            recalc_CV!(fields[i])
            end
        end 

        tmp = zeros(Int64, numinstances)

        for _ in 1:params.Ntherm
            @threads for i in 1:numinstances
                tmp[i] = sweep!(fields[i], rng[i], ϵ, multi_hit)
                ϵ = adjusted_ϵ(ϵ, tmp[i], metro_norm, metro_target_acc)
            end
        end

        if params.starting_Q !== nothing
            @threads for i in 1:numinstances
                instanton!(fields[i], params.starting_Q[i], rng[i], metro_test = false)
            end
        end

        numaccepts = zeros(Int64, numinstances)

        for itrj in 1:params.Nsweeps
            @threads for i in 1:numinstances
                tmp[i] = sweep_meta!(fields[i], cum_bias, rng[i], ϵ, multi_hit)
                ϵ = adjusted_ϵ(ϵ, tmp[i], metro_norm, metro_target_acc)
                numaccepts[i] += tmp[i]
                update_bias!(biases[i], fields[i].CV)
            end

            for i in 1:numinstances
                cum_bias.values .+= biases[i].values 
                clear!(biases[i])
            end

            if params.take_snapshot_every !== nothing
                if itrj%params.take_snapshot_every == 0
                    write_to_file!(cum_bias, open(params.snapshot_dir * "/snapshot_$(itrj).txt", "w"))
                end
            end 

            @threads for i in 1:numinstances
                measurements(itrj, fields[i], meassets[i])
            end
        end

        println_verbose(verbose, "Final ϵ_metro: ", ϵ)

        for i in 1:numinstances
            println_verbose(
                verbose,
                "Acceptance rate ",
                i,
                ": ",
                100 * numaccepts[i] / params.Nsweeps / fields[i].NV / 2 / multi_hit,
                "%"
            )
        end
            
        write_to_file!(cum_bias, cum_bias.fp)

        flush(stdout)
        flush(verbose)
        return nothing
    end

end
