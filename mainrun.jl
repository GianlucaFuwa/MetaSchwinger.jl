module Mainrun
    using Printf 
    using DelimitedFiles
    using InteractiveUtils
    using Dates
    using Base.Threads
    
    import ..System_parameters:Params,Params_set,parameterloading
    import ..Verbose_print:Verbose_,println_verbose,print2file
    import ..Gaugefields:Gaugefield,recalc_S!,recalc_Q!,daction,plaquette
    import ..Metadynamics:Bias_potential,update_bias!,penalty_potential
    import ..Measurements:Measurement_set,measurements,calc_topcharge,calc_weights,build_measurements
    import ..MC:metropolis!,metropolis_meta!,try_swap!,sweep!,sweep_meta!

    import ..System_parameters:physical,meta,sim,mc,meas,system

    function run_sim(filenamein::String)
        filename = filenamein
        include(abspath(filename))
        params_set = Params_set(physical,meta,sim,mc,meas,system)

        run_sim(params_set)

        return nothing
    end

    function run_sim(params_set::Params_set)
        params = parameterloading(params_set)

        if params.meta_enabled
            if params.tempering_enabled
                fields = Vector{Gaugefield}(undef,0)
                biases = Vector{Bias_potential}(undef,0)
                for i=1:params.Ninstances
                push!(fields,Gaugefield(params))
                end
                for i=1:params.Ninstances-1
                push!(biases,Bias_potential(params,i))
                end
                run_temperedsim!(fields,biases,params)
            elseif ~params.parallel_tempering
                field = Gaugefield(params)
                bias = Bias_potential(params,1)
                run_sim!(field,bias,params)
            end
        else
            field = Gaugefield(params)
            run_sim!(field,params)
        end
        return nothing
    end

    function run_sim!(field::Gaugefield,params::Params)
        verbose = Verbose_(params.logfile)
        println_verbose(verbose,"# ",pwd())
        println_verbose(verbose,"# ",Dates.now())
        versioninfo(verbose)

        measset = Measurement_set(params.measure_dir,meas_calls = params.meas_calls)
        ϵ = params.ϵ_metro
        rng = params.randomseeds[1]

        if params.initial == "hot"
            field.g = rand(size(field.g) .- 0.5)*2*2pi
            recalc_S!(field)
            recalc_Q!(field)
        end 

        for therm = 1:params.Ntherm
            sweep!(field,rng,ϵ)
        end
        
        numaccepts = 0

        for trj = 1:params.Nsweeps
            tmp = sweep!(field,rng,ϵ)
            numaccepts += tmp
            # writing logs....#
            if params.veryverbose
            print2file(verbose," ",trj," ",100*numaccepts/trj/field.NV/2,"%"," # trj accrate")
            end
            #-----------------#
            measurements(trj,field,measset)
        end
        println_verbose(verbose,"Acceptance rate: $(100*numaccepts/params.Nsweeps/field.NV/2)%")

        flush(stdout)
        flush(verbose)
        return nothing
    end
    
    # SIMULATION INCLUDING METAD #
    function run_sim!(field::Gaugefield,bias::Bias_potential,params::Params)
        verbose = Verbose_(params.logfile)
        println_verbose(verbose,"# ",pwd())
        println_verbose(verbose,"# ",Dates.now())
        versioninfo(verbose)

        measset = Measurement_set(params.measure_dir,meas_calls = params.meas_calls)
        ϵ = params.ϵ_metro
        rng = params.randomseeds[1]

        if params.initial == "hot"
            field.g = rand(size(field.g) .- 0.5)*2*2pi
            recalc_S!(field)
            recalc_Q!(field)
        end 

        for therm = 1:params.Ntherm
            sweep!(field,rng,ϵ)
        end
        
        numaccepts = 0
        bias_mean = deepcopy(bias.values)
        for trj = 1:params.Nsweeps
            tmp = sweep_meta!(field,bias,rng,ϵ,false)
            numaccepts += tmp
            # writing logs....#
            if params.veryverbose
            print2file(verbose," ",trj," ",100*numaccepts/trj/field.NV/2,"%"," # trj accrate")
            end
            #-----------------#
            bias_mean += bias.values
            seekstart(bias.fp)
            writedlm(bias.fp,[bias.q_vals bias_mean/trj])
            flush(bias)
            measurements(trj,field,measset)
        end
        println_verbose(verbose,"Acceptance rate: $(100*numaccepts/params.Nsweeps/field.NV/2)%")

        bias.values = bias_mean ./ (params.Nsweeps÷2)
        open(params.biasfile,"w") do io
            writedlm(io,[bias.q_vals bias.values])
        end

        q_vals = readdlm(pwd()*"/"*params.measure_dir*"/Continuous_charge.txt",Float64,comments=true)
        weights = calc_weights(q_vals[:,2],bias)
        open(params.weightfile,"w") do io
            writedlm(io,weights)
        end

        flush(stdout)
        flush(verbose)
        return nothing
    end
    
    # SIMULATION INCLUDING METAD & PARALLEL TEMPERING #
    function run_temperedsim!(fields::Vector{Gaugefield},biases::Vector{Bias_potential},params::Params)
        Ninstances = params.Ninstances
        verbose = Verbose_(params.logfile)
        println_verbose(verbose,"# ",pwd())
        println_verbose(verbose,"# ",Dates.now())
        versioninfo(verbose)
        println_verbose(verbose,"# Parallel tempered run with ",params.Ninstances," instances")

        meassets = Vector{Measurement_set}(undef,0)
        for i=1:Ninstances
            push!(meassets,Measurement_set(params.measure_dir,meas_calls=params.meas_calls,instance="_$i"))
        end

        ϵ = params.ϵ_metro
        rng = params.randomseeds

        if params.initial == "hot"
            for i=1:Ninstances
            fields[i].g = rand(size(fields[i]) .- 0.5)*2*2pi
            recalc_Q!(fields[i])
            recalc_S!(fields[i])
            end
        end 

        for therm = 1:params.Ntherm
            @threads for i=1:Ninstances
            sweep!(fields[i],rng[i],ϵ)
            end
        end

        numaccepts = zeros(Int64,Ninstances)
        num_swaps = zeros(Int64,Ninstances-1)
        tmp = zeros(Int64,Ninstances)
        bias_means = []
        for i=1:Ninstances-1
            push!(bias_means,deepcopy(biases[1].values))
        end

        for trj = 1:params.Nsweeps
            tmp[1] = sweep!(fields[1],rng[1],ϵ)
            numaccepts[1] += tmp[1]
            @threads for i=2:Ninstances
                tmp[i] = sweep_meta!(fields[i],biases[i-1],rng[i],ϵ,false)
                numaccepts[i] += tmp[i]
            end
            # writing logs....#
            if params.veryverbose == true
                for i=1:Ninstances
                    print2file(verbose,i," ",trj," ",100*numaccepts[i]/trj/fields[i].NV/2,"%"," # instance trj accrate")
                    print2file(verbose,"------------------------------------------------------")
                end
            end
            #-----------------#

            for i=1:Ninstances-1
                bias_means[i] += biases[i].values
                seekstart(biases[i].fp)
                writedlm(biases[i].fp,[biases[i].q_vals bias_means[i]/trj])
                flush(biases[i])
            end

            if trj%params.swap_every == 0
                for i=Ninstances-1:-1:1
                accept_swap = try_swap!(fields[i],fields[i+1],biases[i],rng[1],false)
                num_swaps[i] += ifelse(accept_swap,1,0)
                if params.veryverbose == true
                    print2file(verbose,i," ⇔ ",i+1," ",100*num_swaps[i]/(trj÷params.swap_every),"%"," # trj swapaccrate")
                end
                end
            end
            print2file(verbose,"------------------------------------------------------")
            @threads for i=1:Ninstances
            measurements(trj,fields[i],meassets[i])
            end
        end # END SWEEPS

        for i=1:Ninstances
        println_verbose(verbose,"Acceptance rate ",i,": ",100*numaccepts[i]/params.Nsweeps/fields[1].NV/2,"%")
        if i<Ninstances
            println_verbose(verbose,"Swap Acceptance rate ",i," ⇔ ",i+1,": ",100*num_swaps[i]/(params.Nsweeps÷params.swap_every),"%")
        end
        end

        for i=1:Ninstances-1
            q_vals = readdlm(pwd()*"/"*params.measure_dir*"/Continuous_charge_$i.txt",Float64,comments=true)
            weights = calc_weights(q_vals[:,2],biases[i])
            open(params.weightfiles[i],"w") do io
                writedlm(io,weights)
            end
            flush(biases[i])
        end

        flush(stdout)
        flush(verbose)
        return nothing
    end

end
