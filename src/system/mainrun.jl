module Mainrun
    using Base.Threads
    using Dates
    using DelimitedFiles
    using InteractiveUtils

    import ..DiracOperators: dirac_operator
    import ..Gaugefields: Gaugefield, random_gaugefield!, recalc_Sg!, recalc_CV!
    import ..Updates: adjusted_ϵ, sweep!, sweep_meta!
    import ..Updates: instanton!, instanton_update!, tempering_swap!
    import ..Metadynamics: BiasPotential, update_bias!, write_to_file!
    import ..Measurements: calc_weights, MeasurementSet, measurements
    import ..Parameters: ParameterSet, ParamSet, parameterloading
    import ..Parameters: physical, meta, param_meta, update, dirac, meas, system
    import ..VerbosePrint: println_verbose, print2file, Verbose_

    function run_sim(filenamein::String)
        filename = filenamein
        include(abspath(filename))
        params_set = ParamSet(physical, meta, param_meta, update, dirac, meas, system)
        run_sim(params_set)
        return nothing
    end

    function run_sim(params_set::ParamSet)
        params = parameterloading(params_set)

        if params.meta_enabled
            if params.tempering_enabled
                fields = Vector{Gaugefield}(undef, 0)
                biases = Vector{BiasPotential}(undef, 0)

                for _ in 1:params.numinstances
                    push!(fields, Gaugefield(params))
                end

                for i in 1:params.numinstances
                    push!(
                        biases,
                        BiasPotential(params, instance = i - 1 + params.no_zero_instance),
                    )
                end

                run_tempered_sim!(fields, biases, params)
            elseif !params.tempering_enabled
                field = Gaugefield(params)
                bias = BiasPotential(params)
                run_sim!(field, bias, params)
            end
        else
            field = Gaugefield(params)
            run_sim!(field, params)
        end

        return nothing
    end

    function run_sim!(field::Gaugefield, params::ParameterSet)
        verbose = Verbose_(params.logfile)
        println_verbose(verbose, "# ", pwd())
        println_verbose(verbose, "# ", Dates.now())
        versioninfo(verbose)

        DiracOperator = dirac_operator(params)

        measset = MeasurementSet(
            params.measure_dir,
            meas_calls = params.meas_calls,
            D = DiracOperator,
        )

        ϵ = params.ϵ_metro
        multi_hit = params.metro_multi_hit
        metro_target_acc = params.metro_target_acc
        metro_norm = 1 / (field.NV * 2 * multi_hit)
        rng = params.randomseeds[1]

        if params.initial == "hot"
            random_gaugefield!(field, rng)
        end

        for _ in 1:params.Ntherm
            tmp = sweep!(field, rng, ϵ, multi_hit)
            ϵ = adjusted_ϵ(ϵ, tmp, metro_norm, metro_target_acc)
        end

        numaccepts = 0
        numaccepts_instanton = 0

        for itrj in 1:params.Nsweeps
            numaccepts += sweep!(field, rng, ϵ, multi_hit)

            if params.instanton_enabled
                Q = ifelse(rand(rng) < 0.5, 1, -1)
                numaccepts_instanton += instanton_update!(field, Q, rng)
            end

            ϵ = adjusted_ϵ(ϵ, tmp, metro_norm, metro_target_acc)

            if params.veryverbose
                println_verbose(
                    verbose,
                    " ",
                    itrj,
                    " ",
                    100 * numaccepts / itrj / field.NV / 2 / multi_hit,
                    "%",
                    " # itrj accrate"
                )
            end

            measurements(itrj, field, measset)
        end

        println_verbose(verbose, "Final ϵ_metro: ", ϵ)
        println_verbose(
            verbose,
            "Acceptance rate: ",
            100 * numaccepts / params.Nsweeps / field.NV / 2 / multi_hit,
            "%"
        )

        if params.instanton_enabled
            println_verbose(
                verbose,
                "Instanton acceptance rate: ",
                100 * numaccepts_instanton / params.Nsweeps / 2,
                "%"
            )
        end

        flush(stdout)
        flush(verbose)
        return nothing
    end

    #=====================================================================================#

    function run_sim!(field::Gaugefield, bias::BiasPotential, params::ParameterSet)
        verbose = Verbose_(params.logfile)
        println_verbose(verbose, "# ", pwd())
        println_verbose(verbose, "# ", Dates.now())
        versioninfo(verbose)

        DiracOperator = dirac_operator(params)

        measset = MeasurementSet(
            params.measure_dir,
            meas_calls = params.meas_calls,
            D = DiracOperator,
        )

        ϵ = params.ϵ_metro
        multi_hit = params.metro_multi_hit
        metro_target_acc = params.metro_target_acc
        metro_norm = 1 / (field.NV * 2 * multi_hit)
        rng = params.randomseeds[1]

        if params.initial == "hot"
            random_gaugefield!(field, rng)
        end

        for _ in 1:params.Ntherm
            tmp = sweep!(field, rng, ϵ, multi_hit)
            ϵ = adjusted_ϵ(ϵ, tmp, metro_norm, metro_target_acc)
        end

        numaccepts = 0
        bias_mean = deepcopy(bias.values)

        for itrj in 1:params.Nsweeps
            numaccepts += sweep_meta!(field, bias, rng, ϵ, multi_hit)
            ϵ = adjusted_ϵ(ϵ, tmp, metro_norm, metro_target_acc)
            update_bias!(bias, field.CV, itrj=itrj)

            if params.veryverbose
                println_verbose(
                    verbose,
                    " ",
                    itrj,
                    " ",
                    100 * numaccepts / itrj / field.NV / 2,
                    "%",
                    " # itrj accrate"
                )
            end

            bias_mean += bias.values
            measurements(itrj, field, measset)
        end

        println_verbose(verbose, "Final ϵ_metro: ", ϵ)
        println_verbose(
            verbose,
            "Acceptance rate: ",
            100 * numaccepts / params.Nsweeps / field.NV / 2 / multi_hit,
            "%"
        )
        println_verbose(verbose, "Exceeded Count: ", bias.exceeded_count)

        copyto!(bias.values, bias_mean / params.Nsweeps)
        write_to_file!(bias)
        flush(bias)

        q_vals = readdlm(
            pwd() * "/" * params.measure_dir * "/Meta_charge.txt",
            Float64,
            comments = true
        )
        weights = calc_weights(q_vals[:,2], bias)

        open(params.weightfiles[1], "w") do io
            writedlm(io, [q_vals[:,1] weights])
        end

        effective_samplesize = round(Int, sum(weights)^2 / sum(weights.^2))
        println_verbose(verbose, "Effective sample size: ", effective_samplesize)

        flush(stdout)
        flush(verbose)
        return nothing
    end

    #=====================================================================================#

    function run_tempered_sim!(
        fields::Vector{Gaugefield},
        biases::Vector{BiasPotential},
        params::ParameterSet
    )
        numinstances = params.numinstances
        verbose = Verbose_(params.logfile)
        println_verbose(verbose, "# ", pwd())
        println_verbose(verbose, "# ", Dates.now())
        versioninfo(verbose)
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

        ϵ = params.ϵ_metro
        multi_hit = params.metro_multi_hit
        metro_target_acc = params.metro_target_acc
        metro_norm = 1 / (fields[1].NV * 2 * multi_hit)
        rng = params.randomseeds

        if params.initial == "hot"
            for i in 1:numinstances
                random_gaugefield!(fields[i], rng[i])
            end
        end

        for _ in 1:params.Ntherm
            @threads for i in 1:numinstances
                numaccepts[i] += sweep!(fields[i], rng[i], ϵ, multi_hit)
                ϵ = adjusted_ϵ(ϵ, numaccepts[i], metro_norm, metro_target_acc)
            end
        end

        if params.starting_Q !== nothing
            @threads for i in 1:numinstances
                instanton!(fields[i], params.starting_Q[i], rng[i], metro_test = false)
            end
        end

        numaccepts = zeros(Int64, numinstances)
        num_swaps = zeros(Int64, numinstances - 1)
        bias_means = []

        for i in 1:numinstances
            push!(bias_means, deepcopy(biases[i].values))
        end

        for itrj in 1:params.Nsweeps

            @threads for i in 1:numinstances
                numaccepts[i] += sweep_meta!(fields[i], biases[i], rng[i], ϵ, multi_hit)
                ϵ = adjusted_ϵ(ϵ, numaccepts[i], metro_norm, metro_target_acc)
                bias_means[i] += biases[i].values
                update_bias!(biases[i], fields[i].CV)
            end

            if itrj%params.swap_every == 0
                for i in numinstances-1:-1:1
                    accept_swap = tempering_swap!(
                        fields[i],
                        fields[i+1],
                        biases[i],
                        biases[i+1],
                        rng[1],
                        actually_swap = false,
                    )
                    num_swaps[i] += accept_swap

                    if params.veryverbose == true && itrj%10000 == 0
                        println_verbose(
                            verbose,
                            num_swaps[i] / (10000 ÷ params.swap_every),
                            " # swapaccrate $i ⇔ $(i+1)",
                        )
                        num_swaps[i] = 0
                    end
                end
            end

            @threads for i in 1:numinstances
                measurements(itrj, fields[i], meassets[i])
            end
        end

        println_verbose(verbose, "# Final ϵ_metro: ", ϵ)

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

        for i in 1+!params.no_zero_instance:numinstances
            println_verbose(verbose, "# Exceeded Count $(i-1): ", biases[i].exceeded_count)
            q_vals = readdlm(
                pwd()*"/$(params.measure_dir)/Meta_charge$(i-!params.no_zero_instance).txt",
                Float64,
                comments = true,
            )
            weights = calc_weights(q_vals[:,2], biases[i])

            open(params.weightfiles[i - !params.no_zero_instance], "w") do io
                writedlm(io, [q_vals[:,1] weights])
            end

            write_to_file!(biases[i])
            effective_samplesize = round(Int, sum(weights)^2 / sum(weights.^2))
            println_verbose(verbose, "# Effective sample size $i: ", effective_samplesize)
        end

        flush(stdout)
        flush(verbose)
        return nothing
    end

end
