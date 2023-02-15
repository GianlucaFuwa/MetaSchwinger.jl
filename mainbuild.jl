module Mainbuild
    using Printf 
    using DelimitedFiles
    using InteractiveUtils
    using Dates
    using Base.Threads
    
    import ..System_parameters:Params,print_parameters,Params_set,parameterloading
    import ..Verbose_print:Verbose_,println_verbose
    import ..Gaugefields:Gaugefield,recalc_S!,recalc_Q!,daction,plaquette
    import ..Metadynamics:Bias_potential,update_bias!,penalty_potential
    import ..Measurements:Measurement_set,calc_topcharge,calc_weights,build_measurements
    import ..MC:metropolis!,metropolis_meta!,try_swap!,sweep!,sweep_meta!

    import ..System_parameters:physical,meta,sim,mc,meas,system

    function run_build(filenamein::String)
        filename = filenamein
        include(abspath(filename))
        params_set = Params_set(physical,meta,sim,mc,meas,system)

        run_build(params_set)

        return nothing
    end

    function run_build(params_set::Params_set)
        params = parameterloading(params_set)
        if params.parallel_tempering
            field_main = Gaugefield(params)
            field_meta = Gaugefield(params)
            bias = Bias_potential(params)
            run_temperedbuild!(field_main,field_meta,bias,params)
        elseif ~params.parallel_tempering
            field = Gaugefield(params)
            bias = Bias_potential(params)
            run_build!(field,bias,params)
        end

        return nothing
    end

    function run_build!(field::Gaugefield,bias::Bias_potential,params::Params)
        verbose = Verbose_(params.logfile)
        println_verbose(verbose,"# ",pwd())
        println_verbose(verbose,"# ",Dates.now())
        versioninfo(verbose)
        println_verbose(verbose,"- - - - - - - - - - - - - - - - - - -")

        measset = Measurement_set(params.measure_dir)
        ϵ = params.ϵ_metro
        rng = params.randomseeds

        if params.initial == "hot"
            field.g = rand(size(field.g) .- 0.5)*2*2pi
        end 

        recalc_S!(field)
        recalc_Q!(field)

        for i = 1:params.Ntherm
            sweep!(field,rng[1],ϵ)
        end

        numaccepts = 0
        for trj = 1:params.Nsweeps
            tmp = sweep_meta!(field,bias,rng[1],ϵ,false)
            numaccepts += tmp
            build_measurements(trj,field,measset)
        end

        open(params.biasfile,"w") do io
            writedlm(io,[bias.q_vals bias.values])
        end    

        println_verbose(verbose,"Acceptance rate: $(100*numaccepts/params.Nsweeps/field.Nx/field.Nt/2)%")
        println_verbose(verbose,"- - - - - - - - - - - - - - - - - - -")

        println_verbose(verbose,"Metapotential has been saved in file \"$(params.biasfile)\"")
        flush(stdout)
        flush(verbose)
        close(verbose.fp)
        return nothing
    end

    function run_temperedbuild!(field_main::Gaugefield,field_meta::Gaugefield,bias::Bias_potential,params::Params)
        verbose = Verbose_(params.logfile)
        println_verbose(verbose,"# ",pwd())
        println_verbose(verbose,"# ",Dates.now())
        versioninfo(verbose)
        println_verbose(verbose,"### This is a parallel tempered run! ###")

        measset_main = Measurement_set(params.measure_dir,meas_calls = params.meas_calls)
        measset_meta = Measurement_set(params.measure_dir_secondary,meas_calls = params.meas_calls)

        ϵ = params.ϵ_metro
        rng = params.randomseeds

        if params.initial == "hot"
            field_main.g = rand(size(field_main.g) .- 0.5)*2*2pi
            field_meta.g = rand(size(field_meta.g) .- 0.5)*2*2pi
        end 

        recalc_S!(field_main)
        recalc_Q!(field_main)
        recalc_S!(field_meta)
        recalc_Q!(field_meta)

        for therm = 1:params.Ntherm
            tmp = @spawn sweep!(field_main,rng[2],ϵ)
            sweep!(field_meta,rng[1],ϵ)
            fetch(tmp)
        end

        numaccepts_main = 0
        numaccepts_meta = 0
        num_swaps = 0

        for trj = 1:params.Nsweeps
            tmp1 = @spawn sweep!(field_main,rng[2],ϵ)
            tmp2 = sweep_meta!(field_meta,bias,rng[1],ϵ,false)
            numaccepts_main += fetch(tmp1)
            numaccepts_meta += tmp2

            if trj%params.swap_every == 0
                accept_swap = try_swap!(field_main,field_meta,bias,rng[1],false)
                num_swaps += ifelse(accept_swap,1,0)
            end
            @sync begin
            @spawn build_measurements(trj,field_main,measset_main)
            build_measurements(trj,field_meta,measset_meta)
            end
        end

        println_verbose(verbose,"Main Acceptance rate: $(100*numaccepts_main/params.Nsweeps/field_main.Nx/field_main.Nt/2)%")
        println_verbose(verbose,"Meta Acceptance rate: $(100*numaccepts_meta/params.Nsweeps/field_meta.Nx/field_meta.Nt/2)%")
        println_verbose(verbose,"Swap Acceptance rate: $(100*num_swaps/(params.Nsweeps÷params.swap_every))%")

        open(params.biasfile,"w") do io
            writedlm(io,[bias.q_vals bias.values])
        end
        println_verbose(verbose,"Metapotential has been saved in file \"$(params.biasfile)\"")

        flush(stdout)
        flush(verbose)
        return nothing
    end

end
