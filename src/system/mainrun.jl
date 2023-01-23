module Mainrun
    using Printf 
    using DelimitedFiles
    using InteractiveUtils
    using Dates
    using Base.Threads:@spawn,nthreads
    
    import ..System_parameters:Params,Params_set,parameterloading
    import ..Verbose_print:Verbose_,println_verbose
    import ..Gaugefields:Gaugefield,recalc_S!,recalc_Q!,instanton,daction,plaquette
    import ..Metadynamics:Bias_potential,update_bias!,penalty_potential
    import ..Measurements:Measurement_set,measurements,calc_topcharge,calc_weights,build_measurements
    import ..MC:metropolis!,metropolis_meta!,instanton_update!,try_swap!,sweep!,sweep_meta!

    import ..System_parameters:physical,meta,sim,system

    function run_sim(filenamein::String)
        filename = filenamein
        include(abspath(filename))
        params_set = Params_set(physical,meta,sim,system)

        run_sim(params_set)

        return nothing
    end

    function run_sim(params_set::Params_set)
        params = parameterloading(params_set)

        if params.parallel_tempering
            @assert nthreads() ≥ 2 "Make sure to enable at least 2 threads when using parallel tempering! Do: 'julia -t 2' or 'julia -t auto'"
            field_main = Gaugefield(params)
            field_meta = Gaugefield(params)
            bias = Bias_potential(params)
            run_temperedsim!(field_main,field_meta,bias,params)
        elseif ~params.parallel_tempering
            field = Gaugefield(params)
            bias = Bias_potential(params)
            run_sim!(field,bias,params)
        end

        return nothing
    end

    function run_sim!(field::Gaugefield,bias::Bias_potential,params::Params)
        verbose = Verbose_(params.logfile)
        println_verbose(verbose,"# ",pwd())
        println_verbose(verbose,"# ",Dates.now())
        versioninfo(verbose)

        measset = Measurement_set(params.measure_dir,meas_calls = params.meas_calls)
        ϵ = params.ϵ
        rng = params.randomseed
        insta1 = instanton(field,1)

        if params.initial == "hot"
            field.g = rand(size(field.g) .- 0.5)*2*2pi
        end 

        if params.meta_runtype == "static"
            static = true
        elseif params.meta_runtype == "dynamic"
            static = false
        end

        recalc_S!(field)
        recalc_Q!(field)
        for trj = 1:params.Ntherm
            for nx = 1:field.Nx
                for nt = 1:field.Nt
                    for d = 1:2
                        metropolis!(field,nx,nt,d,rng,ϵ)
                    end
                end
            end
        end
        if params.weightmode == "from_mean"
            bias_mean = deepcopy(bias.values)
        end
        numaccepts = 0
        for trj = 1:params.Nsweeps
            for eo in [0,1]
                for nx = (eo+1):2:field.Nx
                    for nt = 1:field.Nt
                        accept = metropolis_meta!(field,bias,nx,nt,2,rng,ϵ,static)
                        numaccepts += ifelse(accept,1,0)
                    end
                end
            end
            for eo in [0,1]
                for nt = (eo+1):2:field.Nt
                    for nx = 1:field.Nx
                        accept = metropolis_meta!(field,bias,nx,nt,1,rng,ϵ,static)
                        numaccepts += ifelse(accept,1,0)
                    end
                end
            end
            if trj%params.insta_every == 0
                accept = instanton_update!(field,bias,a,rng,insta1)
                numaccepts += ifelse(accept,1,0)
            end
            measurements(trj,field,measset)
            if params.weightmode == "from_mean"
                bias_mean += bias.values
            end
        end
        if params.weightmode == "from_mean"
            bias.values = bias_mean ./ params.Nsweeps
            open(params.biasfile,"w") do io
                writedlm(io,[bias.q_vals bias.values])
            end
        end
        println_verbose(verbose,"Acceptance rate: $(100*numaccepts/params.Nsweeps/field.Nx/field.Nt/2)%")

        open(params.biasfile,"w") do io
            writedlm(io,[bias.q_vals bias.values])
        end
        println_verbose(verbose,"Metapotential has been saved in file \"$(params.biasfile)\"")
        q_vals = readdlm(pwd()*"/"*params.measure_dir*"/Continuous_charge.txt",Float64,comments=true)
        weights = calc_weights(q_vals[:,2],bias)
        open(params.weightfile,"w") do io
            writedlm(io,weights)
        end
        println_verbose(verbose,"Weights have been saved in file \"$(params.weightfile)\"")
        flush(stdout)
        flush(verbose)
        return nothing
    end

    function run_temperedsim!(field_main::Gaugefield,field_meta::Gaugefield,bias::Bias_potential,params::Params)
        verbose = Verbose_(params.logfile)
        println_verbose(verbose,"# ",pwd())
        println_verbose(verbose,"# ",Dates.now())
        versioninfo(verbose)
        println_verbose(verbose,"### This is a parallel tempered run! ###")

        measset_main = Measurement_set(params.measure_dir,meas_calls = params.meas_calls)
        measset_meta = Measurement_set(params.measure_dir_secondary,meas_calls = params.meas_calls)

        ϵ = params.ϵ
        rng = params.randomseed

        if params.initial == "hot"
            field_main.g = rand(size(field_main.g) .- 0.5)*2*2pi
            field_meta.g = rand(size(field_meta.g) .- 0.5)*2*2pi
        end 

        if params.meta_runtype == "static"
            static = true
        elseif params.meta_runtype == "dynamic"
            static = false
        end

        recalc_S!(field_main)
        recalc_Q!(field_main)
        recalc_S!(field_meta)
        recalc_Q!(field_meta)

        for trj = 1:params.Ntherm
            tmp = @spawn sweep!(field_main,rng,ϵ)
            sweep!(field_meta,rng,ϵ)
            fetch(tmp)
        end

        if params.weightmode == "from_mean"
            bias_mean = deepcopy(bias.values)
        end
        numaccepts_main = 0
        numaccepts_meta = 0
        num_swaps = 0

        for trj = 1:params.Nsweeps
            tmp1 = @spawn sweep!(field_main,rng,ϵ)
            tmp2 = sweep_meta!(field_meta,bias,rng,ϵ,static)
            numaccepts_main += fetch(tmp1)
            numaccepts_meta += tmp2

            if trj%params.swap_every == 0
                accept_swap = try_swap!(field_main,field_meta,bias,rng,static,verbose)
                num_swaps += ifelse(accept_swap,1,0)
            end
            @sync begin
            @spawn measurements(trj,field_main,measset_main)
            measurements(trj,field_meta,measset_meta)
            end
            if params.weightmode == "from_mean"
                bias_mean += bias.values
            end
        end

        if params.weightmode == "from_mean"
            bias.values = bias_mean ./ params.Nsweeps
            open(params.biasfile,"w") do io
                writedlm(io,[bias.q_vals bias.values])
            end
        end
        println_verbose(verbose,"Main Acceptance rate: $(100*numaccepts_main/params.Nsweeps/field_main.Nx/field_main.Nt/2)%")
        println_verbose(verbose,"Meta Acceptance rate: $(100*numaccepts_meta/params.Nsweeps/field_meta.Nx/field_meta.Nt/2)%")
        println_verbose(verbose,"Swap Acceptance rate: $(100*num_swaps/(params.Nsweeps÷params.swap_every))%")

        open(params.biasfile,"w") do io
            writedlm(io,[bias.q_vals bias.values])
        end
        println_verbose(verbose,"Metapotential has been saved in file \"$(params.biasfile)\"")
        q_vals = readdlm(pwd()*"/"*params.measure_dir_secondary*"/Continuous_charge.txt",Float64,comments=true)
        weights = calc_weights(q_vals[:,2],bias)
        open(params.weightfile,"w") do io
            writedlm(io,weights)
        end
        println_verbose(verbose,"Weights have been saved in file \"$(params.weightfile)\"")
        flush(stdout)
        flush(verbose)
        return nothing
    end

end
