module Mainbuild
    using Printf 
    using DelimitedFiles
    using InteractiveUtils
    using Dates
    
    import ..System_parameters:Params,print_parameters,Params_set,parameterloading
    import ..Verbose_print:Verbose_,println_verbose
    import ..Gaugefields:Gaugefield,recalc_S!,recalc_Q!,instanton,daction,plaquette
    import ..Metadynamics:Bias_potential,update_bias!,penalty_potential
    import ..Measurements:Measurement_set,calc_topcharge,calc_weights,build_measurements
    import ..MC:metropolis!,metropolis_meta!,instanton_update!

    import ..System_parameters:physical,meta,sim,system

    function run_build(filenamein::String)
        filename = filenamein
        include(abspath(filename))
        params_set = Params_set(physical,meta,sim,system)

        run_build(params_set)

        return nothing
    end

    function run_build(params_set::Params_set)
        params = parameterloading(params_set)
        field = Gaugefield(params)
        bias = Bias_potential(params)
        run_build!(field,bias,params)

        return nothing
    end

    function run_build!(field::Gaugefield,bias::Bias_potential,params::Params)
        verbose = Verbose_(params.logfile)
        println_verbose(verbose,"# ",pwd())
        println_verbose(verbose,"# ",Dates.now())
        versioninfo(verbose)
        println_verbose(verbose,"- - - - - - - - - - - - - - - - - - -")

        measset = Measurement_set(params.measure_dir)
        ϵ = params.ϵ
        rng = params.randomseed
        
        insta1 = instanton(field,1)

        if params.initial == "hot"
            field.g = rand(size(field.g) .- 0.5)*2*2pi
        end 

        static = false

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
            build_measurements(trj,field,measset)
        end
        open(params.biasfile,"w") do io
            writedlm(io,[bias.q_vals bias.values])
        end
        println_verbose(verbose,"Acceptance rate: $(100*numaccepts/params.Nsweeps/field.Nx/field.Nt/2)%")
        println_verbose(verbose,"- - - - - - - - - - - - - - - - - - -")

        open(params.biasfile,"w") do io
            writedlm(io,[bias.q_vals bias.values])
        end
        println_verbose(verbose,"Metapotential has been saved in file \"$(params.biasfile)\"")
        flush(stdout)
        flush(verbose)
        return nothing
    end

end
