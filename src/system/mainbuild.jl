module Mainbuild
    using Base.Threads
    using Dates
    using DelimitedFiles
    using InteractiveUtils

    import ..BiasModule: Metadynamics, OPES, clear!, parametric_to_bias!, update_bias!,
        dump_state_to_file, MSE
    import ..DiracOperators: dirac_operator
    import ..Gaugefields: Gaugefield, recalc_Sg!, recalc_CV!
    import ..Updates: Updatemethod, update!
    import ..Updates: instanton!, instanton_update!, tempering_swap!
    import ..Measurements: MeasurementSet, calc_weights, measurements
    import ..Parameters: ParameterSet, ParamSet, parameterloading
    import ..Parameters: physical, meta, param_meta, opes, update, dirac, meas, system
    import ..VerbosePrint: print2file, println_verbose, Verbose_

    function run_build(filenamein::String)
        filename = filenamein
        include(abspath(filename))
        params_set = ParamSet(physical, meta, param_meta, opes, update, dirac, meas, system)
        run_build(params_set)
        return nothing
    end

    function run_build(params_set::ParamSet)
        params = parameterloading(params_set)
        verbose = Verbose_(params.logfile)
        println_verbose(verbose, "# ", pwd())
        println_verbose(verbose, "# ", Dates.now())
        versioninfo(verbose)

        if params.tempering_enabled
            numinstances = params.numinstances
            field1 = Gaugefield(params)
            updatemethod1 = Updatemethod(params, field1)

            fields = Vector{Gaugefield}(undef, numinstances)
            updatemethod = Vector{typeof(updatemethod1)}(undef, numinstances)

            for i in 1:numinstances
                fields[i] = deepcopy(field1)
                updatemethod[i] = deepcopy(updatemethod1)
            end

            cum_bias = params.opes_enabled ? OPES(params) : Metadynamics(params)
            run_tempered_build!(fields, cum_bias, updatemethod, params)
        elseif !params.tempering_enabled
            field = Gaugefield(params)
            bias = params.opes_enabled ? OPES(params) : Metadynamics(params)
            updatemethod = Updatemethod(params, field)
            run_build!(field, bias, updatemethod, verbose, params)
        end

        return nothing
    end

    function run_build!(field, bias, updatemethod, verbose, params)
        measset = MeasurementSet(params.measure_dir, meas_calls = params.meas_calls)
        rng = params.randomseeds[1]
        Nsweeps = params.Nsweeps

        if params.initial == "hot"
            field.U = rand(size(field.U) .- 0.5) * 2 * 2π
            recalc_Sg!(field)
            recalc_CV!(field)
        end

        for _ in 1:params.Ntherm
            update!(updatemethod, field, rng)
        end

        if typeof(bias) == Metadynamics
            bias.parametric && parametric_to_bias!(bias)
        end

        numaccepts = 0

        for itrj in 1:Nsweeps
            numaccepts += update!(updatemethod, field, rng; bias=bias)
            update_bias!(bias, field.CV; itrj=itrj)

            if params.veryverbose
                println_verbose(verbose, itrj, " ", 100*numaccepts/itrj, "% # trj accrate")
            end

            measurements(itrj, field, measset)
            # MSEᵢ, STDᵢ = MSE(bias)
            # itrj%100==0 && println_verbose(verbose, MSEᵢ, "\t", STDᵢ, "\t# $itrj")
            # if MSEᵢ < 2.5
            #     println_verbose(verbose, ">>Bias has converged after $(itrj) iterations\n")
            #     break
            # end
            flush(verbose)
        end

        dump_state_to_file(bias, "final")

        println_verbose(verbose, "# Acceptance rate: ", 100*numaccepts/Nsweeps, "%")
        flush(stdout)
        flush(verbose)
        return nothing
    end

    function run_tempered_build!(fields::Vector{Gaugefield}, cum_bias, updatemethod, params)
        numinstances = params.numinstances
        verbose = Verbose_(params.logfile)
        println_verbose(verbose, "# ", pwd())
        println_verbose(verbose, "# ", Dates.now())
        versioninfo(verbose)
        println_verbose(verbose, "Parallel tempered run with $numinstances MetaD instances")

        meassets = Vector{MeasurementSet}(undef, 0)

        for i in 1:numinstances
            push!(
                meassets,
                MeasurementSet(params.measure_dir, meas_calls=params.meas_calls, instance="_$i")
            )
        end

        rng = params.randomseeds

        if params.initial == "hot"
            for i in 1:numinstances
            fields[i].U = rand(size(fields[i].U) .- 0.5) * 2 * 2π
            recalc_Sg!(fields[i])
            recalc_CV!(fields[i])
            end
        end

        if params.starting_Q !== nothing
            println_verbose(verbose, ">> Starting at Qs: ")
            for i in 1:numinstances
                instanton!(fields[i], params.starting_Q[i], rng[i]; metro_test=false)
                println_verbose(verbose, "$(params.starting_Q[i])")
            end
        end

        for _ in 1:params.Ntherm
            for i in 1:numinstances
                update!(updatemethod[i], fields[i], rng[i])
            end
        end

        if typeof(cum_bias) == Metadynamics
            cum_bias.parametric && parametric_to_bias!(cum_bias)
        end

        numaccepts = zeros(Float64, numinstances)
        CVs = zeros(Float64, numinstances)

        for itrj in 1:params.Nsweeps
            for i in 1:numinstances
                numaccepts[i] += update!(updatemethod[i], fields[i], rng[i]; bias=cum_bias)
            end

            for i in 1:numinstances
                CVs[i] = fields[i].CV
            end

            update_bias!(cum_bias, CVs; itrj=itrj)

            if params.take_snapshot_every !== nothing
                if itrj%params.take_snapshot_every == 0
                    write_to_file(cum_bias, open(params.snapshot_dir * "/snapshot_$(itrj).txt", "w"))
                end
            end

            for i in 1:numinstances
                measurements(itrj, fields[i], meassets[i])
            end

            if itrj%1000 == 0
                MSEᵢ, STDᵢ = MSE(cum_bias)
                println_verbose(verbose, MSEᵢ, "\t", STDᵢ, "\t# $itrj")
            end
            flush(verbose)
        end

        dump_state_to_file(cum_bias, "")

        for i in 1:numinstances
            println_verbose(verbose, "# Acc. rate $i: $(100*numaccepts[i]/params.Nsweeps)%")
        end

        flush(stdout)
        flush(verbose)
        return nothing
    end

end
