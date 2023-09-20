module Mainrun
    using Base.Threads
    using Dates
    using DelimitedFiles
    using InteractiveUtils

    import ..BiasModule: Metadynamics, OPES, update_bias!
    import ..DiracOperators: dirac_operator
    import ..Gaugefields: Gaugefield, random_gaugefield!, recalc_Sg!, recalc_CV!
    import ..Updates: Updatemethod, update!
    import ..Updates: instanton!, instanton_update!, tempering_heatbath!, tempering_swap!
    import ..Measurements: MeasurementSet, calc_weights, measurements
    import ..Parameters: ParameterSet, ParamSet, parameterloading
    import ..Parameters: physical, meta, param_meta, opes, update, dirac, meas, system
    import ..VerbosePrint: println_verbose, print2file, Verbose_

    function run_sim(filenamein::String)
        filename = filenamein
        include(abspath(filename))
        params_set = ParamSet(physical, meta, param_meta, opes, update, dirac, meas, system)
        run_sim(params_set)
        return nothing
    end

    function run_sim(params_set::ParamSet)
        params = parameterloading(params_set)
        verbose = Verbose_(params.logfile)
        println_verbose(verbose, "# ", pwd())
        println_verbose(verbose, "# ", Dates.now())
        versioninfo(verbose)
        println_verbose(verbose)

        dirac = dirac_operator(params)
        if params.meta_enabled
            if params.tempering_enabled
                fields = Vector{Gaugefield}(undef, 0)
                biases = Vector{Metadynamics}(undef, 0)

                push!(fields, Gaugefield(params; instance=params.no_zero_instance))
                updatemethod₁ = Updatemethod(params, fields[1])
                updatemethod = Vector{typeof(updatemethod₁)}(undef, 0)
                push!(updatemethod, updatemethod₁)

                for _ in 2:params.numinstances
                    deepcopy(updatemethod₁)
                end

                for i in 1:params.numinstances
                    push!(
                        biases,
                        Metadynamics(params; instance=i-1+params.no_zero_instance),
                    )
                end

                run_tempered_sim!(fields, biases, updatemethod, verbose, params)
            elseif !params.tempering_enabled
                field = Gaugefield(params)
                bias = Metadynamics(params)
                updatemethod = Updatemethod(params, field)
                run_sim!(field, bias, updatemethod, dirac, verbose, params)
            end
        elseif params.opes_enabled
            field = Gaugefield(params)
            bias = OPES(params)
            updatemethod = Updatemethod(params, field)
            run_sim!(field, bias, updatemethod, dirac, verbose, params)
        else
            field = Gaugefield(params)
            updatemethod = Updatemethod(params, field)
            run_sim!(field, updatemethod, dirac, verbose, params)
        end

        return nothing
    end

    function run_sim!(field, updatemethod, dirac, verbose, params)
        measset = MeasurementSet(
            params.measure_dir,
            meas_calls = params.meas_calls,
            D = dirac,
        )

        rng = params.randomseeds[1]

        if params.initial == "hot"
            random_gaugefield!(field, rng)
        end

        for _ in 1:params.Ntherm
            update!(updatemethod, field, rng)
        end

        numaccepts = 0
        numaccepts_instanton = 0

        for itrj in 1:params.Nsweeps
            numaccepts += update!(updatemethod, field, rng)

            if params.instanton_enabled
                Q = ifelse(rand(rng)<0.5, 1, -1)
                numaccepts_instanton += instanton_update!(field, Q, rng)
            end

            measurements(itrj, field, measset)
        end

        println_verbose(
            verbose,
            "Acceptance rate: ",
            100 * numaccepts / params.Nsweeps,
            "%"
        )

        if params.instanton_enabled
            println_verbose(
                verbose,
                "Instanton acceptance rate: ",
                100 * numaccepts_instanton / params.Nsweeps,
                "%"
            )
        end

        flush(stdout)
        flush(verbose)
        return nothing
    end

    #=====================================================================================#

    function run_sim!(field, bias, updatemethod, dirac, verbose, params)
        measset = MeasurementSet(
            params.measure_dir,
            meas_calls = params.meas_calls,
            D = dirac,
        )

        rng = params.randomseeds[1]

        if params.initial == "hot"
            random_gaugefield!(field, rng)
        end

        for _ in 1:params.Ntherm
            update!(updatemethod, field, rng)
        end

        numaccepts = 0
        if typeof(bias) <: Metadynamics
            bias_mean = deepcopy(bias.values)
        end

        for itrj in 1:params.Nsweeps
            numaccepts += update!(updatemethod, field, rng; bias=bias)
            update_bias!(bias, field.CV, itrj=itrj)

            if params.veryverbose
                println_verbose(
                    verbose,
                    " ",
                    itrj,
                    " ",
                    100 * numaccepts / itrj,
                    "%",
                    " # itrj accrate"
                )
            end

            typeof(bias) <: Metadynamics && (bias_mean += bias.values)
            measurements(itrj, field, measset)
        end

        println_verbose(
            verbose,
            "Acceptance rate: ",
            100 * numaccepts / params.Nsweeps,
            "%"
        )

        if typeof(bias) <: Metadynamics
            copyto!(bias.values, bias_mean / params.Nsweeps)

            q_vals, _ = readdlm(
                pwd() * "/" * params.measure_dir * "/Meta_charge.txt",
                Float64,
                comments = true,
                header = true,
            )
            weights = calc_weights(q_vals[:,2], bias)

            open(params.weightfiles[1], "w") do io
                writedlm(io, [q_vals[:,1] weights])
            end

            effective_samplesize = round(Int, sum(weights)^2 / sum(weights.^2))
            println_verbose(verbose, "Effective sample size: ", effective_samplesize)
        end

        flush(stdout)
        flush(verbose)
        return nothing
    end

    #=====================================================================================#

    function run_tempered_sim!(
        fields::Vector{Gaugefield},
        biases::Vector{Metadynamics},
        updatemethod,
        verbose,
        params,
    )
        numinstances = params.numinstances
        println_verbose(
            verbose,
            "# Parallel tempered run with ",
            Int(!params.no_zero_instance),
            " regular instance and ",
            numinstances - !params.no_zero_instance,
            " MetaD instances"
        )

        meassets = Vector{MeasurementSet}(undef, 0)

        for i in 1:numinstances
            push!(
                meassets,
                MeasurementSet(params.measure_dir, meas_calls = params.meas_calls,
                instance = "_$(i - !params.no_zero_instance)")
            )
        end

        rng = params.randomseeds
        tempering_heatbath = false

        if params.initial == "hot"
            for i in 1:numinstances
                random_gaugefield!(fields[i], rng[i])
            end
        end

        for _ in 1:params.Ntherm
            for i in 1:numinstances
                update!(updatemethod[i], fields[i], rng[i])
            end
        end

        if params.starting_Q !== nothing
            for i in 1:numinstances
                instanton!(fields[i], params.starting_Q[i], rng[i], metro_test=false)
            end
        end

        numaccepts = zeros(Int64, 8numinstances)
        num_swaps = zeros(Int64, numinstances - 1)
        bias_means = []

        for i in 1:numinstances
            push!(bias_means, deepcopy(biases[i].values))
        end

        for itrj in 1:params.Nsweeps

            for i in 1:numinstances
                numaccepts[i] = update!(updatemethod[i], fields[i], rng[i]; bias=biases[i])
                bias_means[i] += biases[i].values
                update_bias!(biases[i], fields[i].CV)
            end

            if itrj%params.swap_every == 0
                if tempering_heatbath#params.ptmetad_heatbath
                    num_swaps[1] += tempering_heatbath!(fields, biases, rng[1])
                    if itrj%5000 == 0
                        println_verbose(
                            verbose,
                            num_swaps[1] / (5000÷params.swap_every),
                            " # swapaccrate",
                        )
                        num_swaps[1] = 0
                    end
                else
                    for i in numinstances-1:-1:1
                        num_swaps[i] += tempering_swap!(
                            fields[i],
                            fields[i+1],
                            biases[i],
                            biases[i+1],
                            rng[1],
                        )
                        if itrj%5000 == 0
                            println_verbose(
                                verbose,
                                num_swaps[i] / (5000÷params.swap_every),
                                " # swapaccrate $i ⇔ $(i+1)",
                            )
                            num_swaps[i] = 0
                        end # END PRINT
                    end # END SWAP LOOP
                end # END IF HEATBATH
            end # END IF SWAP

            for i in 1:numinstances
                measurements(itrj, fields[i], meassets[i])
            end
        end

        if tempering_heatbath
            # println_verbose(verbose, "# Swap Acceptance rate: $(num_swaps[1]/)")
        else
            for i in 1:numinstances
                println_verbose(
                    verbose,
                    "# Acceptance rate ",
                    i,
                    ": ",
                    100 * numaccepts[i] / params.Nsweeps / fields[i].NV / 2 / multi_hit,
                    "%"
                )
                if i < numinstances
                    println_verbose(
                        verbose,
                        "# Swap Acceptance rate ",
                        i,
                        " ⇔ ",
                        i + 1,
                        ": ",
                        100 * num_swaps[i] / (params.Nsweeps ÷ params.swap_every),
                        "%"
                    )
                end
            end
        end

        for i in 1+!params.no_zero_instance:numinstances
            q_vals, _ = readdlm(
                pwd()*"/$(params.measure_dir)/Meta_charge_$(i-!params.no_zero_instance).txt",
                header = true,
            )
            weights = calc_weights(q_vals[:,2], biases[i])

            open(params.weightfiles[i - !params.no_zero_instance], "w") do io
                writedlm(io, [q_vals[:,1] weights])
            end

            effective_samplesize = round(Int, sum(weights)^2 / sum(weights.^2))
            println_verbose(verbose, "# Effective sample size $i: ", effective_samplesize)
        end

        flush(stdout)
        flush(verbose)
        return nothing
    end

end
