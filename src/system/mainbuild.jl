module Mainbuild
    using DelimitedFiles
    using InteractiveUtils
    using Dates
    
    import ..System_parameters: Params,Params_set,parameterloading
    import ..Verbose_print: Verbose_,println_verbose,print2file
    import ..Gaugefields: Gaugefield,recalc_Sg!,recalc_CV!
    import ..Metadynamics: Bias_potential,update_bias!,write_to_file!,parametric_to_bias!
    import ..Measurements: Measurement_set,measurements,calc_weights
    import ..Local: sweep!,sweep_meta!
    import ..Tempering: tempering_swap!

    import ..System_parameters:physical,meta,param_meta,sim,mc,meas,system

    function run_build(filenamein::String)
        filename = filenamein
        include(abspath(filename))
        params_set = Params_set(physical,meta,param_meta,sim,mc,meas,system)

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

        measset = Measurement_set(params.measure_dir,meas_calls = params.meas_calls)
        ϵ = params.ϵ_metro
        rng = params.randomseeds[1]
        Nsweeps = params.Nsweeps

        if params.initial == "hot"
            field.U = rand(size(field.U) .- 0.5)*2*2pi
            recalc_Sg!(field)
            recalc_CV!(field)
        end 

        for therm = 1:params.Ntherm
            sweep!(field,rng,ϵ)
        end
        #
        if bias.parametric == true
            parametric_to_bias!(bias)
        end
        #
        numaccepts = 0
        for itrj = 1:Nsweeps
            tmp = sweep_meta!(field,bias,rng,ϵ)
            numaccepts += tmp
            update_bias!(bias, field.CV, itrj=itrj)
            # writing logs....#
            if params.veryverbose
            println_verbose(verbose," ",itrj," ",100*numaccepts/trj/field.NV/2,"%"," # trj accrate")
            end
            #-----------------#
            measurements(itrj,field,measset)
        end
        println_verbose(verbose,"Acceptance rate: ",100*numaccepts/Nsweeps/field.NV/2,"%")
        println_verbose(verbose,"Exceeded Count: ",bias.exceeded_count)

        write_to_file!(bias)
        flush(bias)

        flush(stdout)
        flush(verbose)
        return nothing
    end

end
