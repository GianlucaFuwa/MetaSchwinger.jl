module System_parameters
	using Random
    using DelimitedFiles

    export Params
	
	printlist_physical = [
        "N",
        "β",
    ]

	printlist_meta = [
        "meta_enabled",
        "parametric",
        "symmetric",
        "well_tempered",
        "CVlims",
        "bin_width",
        "w",
        "k",
        "is_static",
        "ΔT",
        "take_snapshot_every",
    ]

    printlist_param_meta = [
        "parametric",
        "potential_parameters",
        "lower_bounds",
        "upper_bounds",
        "batchsize",
        "testfun",
        "minimizer",
    ]

	printlist_sim = [
        "Ntherm",
        "Nsweeps",
        "initial",
        "starting_Q",
        "instanton_enabled",
        "tempering_enabled",
        "Ninstances",
        "no_zero_instance",
        "swap_every",
    ]

    printlist_mc = [
        "update_method",
        "ϵ_metro",
        "multi_hit",
        "metro_target_acc",
        "ϵ_hmc",
        "hmc_steps",
    ]

    printlist_dirac = [
        "Dirac_operator",
        "mass",
        "BC",
        ]

    printlist_meas = [
        "meas_calls",
    ]

    printlist_system = [
        "veryverbose",
        "randomseeds",
        "logdir",
        "logfile",
        "loadfile",
        "measure_basedir",
        "measure_dir",
        "savebias_dir",
        "biasfile",
        "usebiases",
        "weightfile",
        "snapshot_dir",
    ]

    const printlists = [
        printlist_physical,
        printlist_meta,
        printlist_param_meta,
        printlist_sim,
        printlist_mc,
        printlist_dirac,
        printlist_meas,
        printlist_system,
    ]

    const printlists_header = [
        "# Physical Settings ",
        "# Metadynamics Settings",
        "# Parametric Settings",
        "# Simulation Settings",
        "# MC Settings",
        "# Dirac Settings",
        "# Measurement Settings",
        "# System Settings",
    ]

	defaultmeasures = Array{Dict,1}(undef, 2)

	for i in 1:length(defaultmeasures)
		defaultmeasures[i] = Dict()
	end

	defaultmeasures[1]["methodname"] = "Meta_charge"
    defaultmeasures[1]["measure_every"] = 1
    defaultmeasures[2]["methodname"] = "Topological_charge"
    defaultmeasures[2]["measure_every"] = 1

	physical = Dict()
	meta = Dict()
    param_meta = Dict()
	sim = Dict()
    mc = Dict()
    dirac = Dict()
    meas = Dict()
	system = Dict()

    meta["meta_enabled"] = true
    meta["parametric"] = false
    meta["is_static"] = [false]
    meta["well_tempered"] = false
    meta["take_snapshot_every"] = nothing

    param_meta["lower_bounds"] = nothing
    param_meta["upper_bounds"] = nothing
    param_meta["batchsize"] = nothing
    param_meta["testfun"] = nothing
    param_meta["minimizer"] = nothing

    sim["starting_Q"] = nothing
    sim["instanton_enabled"] = false
    sim["tempering_enabled"] = false
    sim["no_zero_instance"] = false

    mc["update_method"] = "Local"
    mc["multi_hit"] = 1
    mc["metro_target_acc"] = 0.7

    dirac["Dirac_operator"] = nothing
    dirac["mass"] = nothing
    dirac["BC"] = nothing

	meas["meas_calls"] = defaultmeasures

    system["veryverbose"] = false
    system["randomseeds"] = [Random.Xoshiro()]
	
	mutable struct ParamSet
	physical::Dict
	meta::Dict
    param_meta::Dict
	sim::Dict
    mc::Dict
    dirac::Dict
    meas::Dict
	system::Dict
	
        function ParamSet(physical, meta, param_meta, sim, mc, dirac, meas, system)
            return new(physical, meta, param_meta, sim, mc, dirac, meas, system)
        end
    end
	
    function make_parameters(physical, meta, param_meta, sim, mc, dirac, meas, system)
        return ParamSet(physical, meta, param_meta, sim, mc, dirac, meas, system)
    end
    
	struct Params
		N::Tuple
		β::Float64

        meta_enabled::Bool
        parametric::Union{Nothing, Bool}
        symmetric::Union{Nothing, Bool}
		CVlims::Union{Nothing, NTuple{2,Float64}}
		bin_width::Union{Nothing, Float64}
		w::Union{Nothing, Float64}
		k::Union{Nothing, Float64}
        is_static::Union{Nothing, Vector{Bool}}
        well_tempered::Union{Nothing, Bool}
        ΔT::Union{Nothing, Float64}
        take_snapshot_every::Union{Nothing, Int64}

        potential_parameters::Union{Nothing, Vector{Float64}}
        lower_bounds::Union{Nothing, Vector{Float64}}
        upper_bounds::Union{Nothing, Vector{Float64}}
        batchsize::Union{Nothing, Int64}
        testfun::Union{Nothing, String}
        minimizer::Union{Nothing, String}

		Ntherm::Int64
		Nsweeps::Int64
		initial::String
        starting_Q::Union{Nothing, Vector{Int64}}
        instanton_enabled::Bool
        tempering_enabled::Union{Nothing, Bool}
        Ninstances::Union{Nothing, Int64}
        no_zero_instance::Union{Nothing, Bool}
        swap_every::Union{Nothing, Int64}

        update_method::String
		ϵ_metro::Union{Float64, Nothing}
        multi_hit::Union{Float64, Nothing}
        metro_target_acc::Union{Float64, Nothing}

        Dirac_operator::Union{String, Nothing}
        mass::Union{Float64, Nothing}
        BC::Union{Vector{Int64}, Nothing}

        meas_calls::Array{Dict,1}

        veryverbose::Bool
		randomseeds::Vector{Xoshiro}
		logdir::String
		logfile::String
		loadfile::IOStream
		measure_dir::String
		savebias_dir::Union{Nothing, String}
		biasfiles::Union{Nothing, Vector{String}}
		usebiases::Union{Nothing, Vector{Union{Nothing, String}}}
        weightfiles::Union{Nothing, Vector{String}}
        snapshot_dir::Union{Nothing, String}

		function Params(physical, meta, param_meta, sim, mc, dirac, meas, system)
			N = physical["N"]
			β = physical["β"]

            measure_dir = system["measure_basedir"] * "/" * system["measure_dir"] * "/raw_data"

            meta_enabled = meta["meta_enabled"]

            if meta_enabled == true
                parametric = meta["parametric"]
                symmetric = meta["symmetric"]
                CVlims = meta["CVlims"]
                bin_width = meta["bin_width"]
                w = meta["w"]
                k = meta["k"]
                is_static = meta["is_static"]
                well_tempered = meta["well_tempered"]
                ΔT = well_tempered==true ? meta["ΔT"] : nothing
                take_snapshot_every = meta["take_snapshot_every"]
                potential_parameters = parametric==true ? param_meta["potential_parameters"] : nothing
                lower_bounds = parametric==true ? param_meta["lower_bounds"] : nothing
                upper_bounds = parametric==true ? param_meta["upper_bounds"] : nothing
                batchsize = parametric==true ? param_meta["batchsize"] : nothing
                testfun = parametric==true ? param_meta["testfun"] : nothing
                minimizer = parametric==true ? param_meta["minimizer"] : nothing

                savebias_dir = system["savebias_dir"]
                
                biasfiles = []
                usebiases = []
                weightfiles = []
                tempering_enabled = sim["tempering_enabled"]

                if tempering_enabled == true
                    swap_every = sim["swap_every"]
                    Ninstances = sim["Ninstances"]
                    no_zero_instance = sim["no_zero_instance"]

                    if haskey(system,"usebiases")
                        usebiases = system["usebiases"]
                    else   
                        push!(usebiases,nothing)
                    end

                    for i in 1:Ninstances - !no_zero_instance
                        push!(biasfiles, pwd() * "/" * savebias_dir * "/" * system["biasfile"]*"_$i.txt")
                        push!(weightfiles, pwd() * "/" * measure_dir * "/Weights_$i.txt")
                    end
                else # IF NO TEMPERING
                    push!(biasfiles, pwd() * "/" * savebias_dir * "/" * system["biasfile"] * ".txt")
                    push!(weightfiles, pwd() * "/" * measure_dir * "/Weights.txt")
                    
                    if haskey(system, "usebiases")
                        push!(usebiases, system["usebiases"][1])
                    else   
                        push!(usebiases, nothing)
                    end
                end # END IF TEMPERING

                if isdir(savebias_dir) == false
                    mkpath(savebias_dir)
                end

                if take_snapshot_every !== nothing
                    snapshot_dir = savebias_dir * "/" * system["snapshot_dir"]

                    if isdir(snapshot_dir) == false
                        mkpath(snapshot_dir)
                    end
                end
            else # IF NO META
                parametric = nothing
                potential_parameters = nothing
                lower_bounds = nothing
                upper_bounds = nothing
                batchsize = nothing
                testfun = nothing
                minimizer = nothing
                symmetric = nothing
                CVlims = nothing
                bin_width = nothing
                w = nothing
                k = nothing
                is_static = nothing
                well_tempered = nothing
                ΔT = nothing
                take_snapshot_every = nothing
                savebias_dir = nothing
                biasfiles = nothing
                usebiases = nothing
                weightfiles = nothing
                snapshot_dir = nothing
                tempering_enabled = nothing
            end # END IF META

            starting_Q = sim["starting_Q"]
            instanton_enabled = sim["instanton_enabled"]
            Ntherm = sim["Ntherm"]
			Nsweeps = sim["Nsweeps"]
            initial = sim["initial"]

            update_method = mc["update_method"]

            if update_method == "Local" || update_method == "Local-Meta"
			    ϵ_metro = mc["ϵ_metro"]
                multi_hit = mc["multi_hit"]
                metro_target_acc = mc["metro_target_acc"]
            else
                error("Only local update method supported")
            end

            Dirac_operator = dirac["Dirac_operator"]
            mass = dirac["mass"]
            BC = dirac["BC"]

            meas_calls = meas["meas_calls"]

            if tempering_enabled == true
                Ninstances = sim["Ninstances"]
                no_zero_instance = sim["no_zero_instance"]
                swap_every = sim["swap_every"]
            else
                Ninstances = nothing
                no_zero_instance = nothing
                swap_every = nothing
            end

            veryverbose = system["veryverbose"]
			randomseeds = system["randomseeds"]
			logdir = system["logdir"]

			if isdir(logdir) == false
				mkpath(logdir)
			end

			logfile = pwd() * "/" * logdir * "/" * system["logfile"] * ".txt"
			loadfile = open(logfile, "a")

            if isdir(measure_dir) == false
                mkpath(measure_dir)
            end
            
			return new(
				N, β,
                meta_enabled, parametric, symmetric, CVlims, bin_width, w, k, is_static, well_tempered, ΔT, take_snapshot_every,
                potential_parameters, lower_bounds, upper_bounds, batchsize, testfun, minimizer,
                Ntherm, Nsweeps, initial, starting_Q, instanton_enabled, tempering_enabled, Ninstances, no_zero_instance, swap_every,
                update_method, ϵ_metro, multi_hit, metro_target_acc,
                Dirac_operator, mass, BC,
                meas_calls,
                veryverbose, randomseeds, logdir, logfile, loadfile, measure_dir, savebias_dir, biasfiles, usebiases, weightfiles, snapshot_dir,
            )
		end

		function Params(params_set::ParamSet)
			return Params(
                params_set.physical,
                params_set.meta,
                params_set.param_meta,
                params_set.sim,
                params_set.mc,
                params_set.dirac,
                params_set.meas,
                params_set.system,
            )
		end
	end
    
    function set_params(dict, string)
        if haskey(dict, string)
            return dict[string]
        else
            error("No parameters given. \n")
        end
    end
	
	function make_parametersdict(p::T) where T<:Any
		pnames = fieldnames(T)
		pdict = Dict()

		for i in eachindex(pnames)
			pdict[String(pnames[i])] = getfield(p,pnames[i])
		end

		return pdict, pnames
	end

	function print_parameters_file(p)
		filename = p.logfile*"_parameters.jl"
		fp = open(filename, "w")
		println(fp, "# - - parameters - - - - - - - - ")
		pdict, = make_parametersdict(p)

		for param in pdict
			if typeof(param[2]) == String
				println(fp, "$(param[1]) = \"$(param[2])\"")
			else
				println(fp, "$(param[1]) = $(param[2])")
			end
		end

		println(fp, "# - - - - - - - - - - - - - - - -")
		close(fp)
        return nothing
	end
	
    function setprint(fp, fp2, string)
        println(fp, string)

        if fp2 !== nothing
            println(fp2, string)
        end

        println(string)
        return nothing
    end
    
    function get_stringfromkey(key)
        if typeof(key) == String
            string = "\"$key\""
        else
            string = "$key"
		end

        return string
	end

	function get_header(params_set::ParamSet, inputname)
		if haskey(params_set.system, inputname)
            return "system", params_set.system[inputname]
        elseif haskey(params_set.physical, inputname)
            return "physical", params_set.physical[inputname]
        elseif haskey(params_set.meta, inputname)
            return "meta", params_set.meta[inputname]
        elseif haskey(params_set.param_meta, inputname)
            return "param_meta", params_set.param_meta[inputname]
        elseif haskey(params_set.mc, inputname)
            return "mc", params_set.mc[inputname]
        elseif haskey(params_set.sim, inputname)
            return "sim", params_set.sim[inputname]
        elseif haskey(params_set.dirac, inputname)
            return "dirac", params_set.dirac[inputname]
        elseif haskey(params_set.meas, inputname)
            return "meas", params_set.meas[inputname]
        else
            return nothing, nothing
        end
    end

	function print_measurementinfo(fp, fp2, key)
        string = "measurement[\"measurement_methods\"] = Dict[ "
        setprint(fp, fp2, string)

        for (i, data) in enumerate(key)
            string = "\t Dict{Any,Any}(\"methodname\" => \"$(data["methodname"])\","
            setprint(fp, fp2, string)
            count = 0
            for (name, key_i) in data
                if name != "methodname"
                    count += 1
                    paramstring = get_stringfromkey(key_i)
                    string = "    \"$(name)\" => " * paramstring

                    if count != length(data) - 1
                        string *= ","
                    end

                    setprint(fp, fp2, string)
                end
            end

            string = "  )"

            if i != length(key)
                string *= ","
            end 

            setprint(fp, fp2, string)
        end

        string = "]"
        setprint(fp, fp2, string)
        return nothing
    end

    function print_parameters_list(params_set::ParamSet, p = nothing; filename = nothing)
        if filename === nothing 
            @assert p !== nothing "wrong input!"

            filename = p.logfile * "_parameters.jl"
            fp2 = p.loadfile
        else
            fp2 = nothing
        end

        fp = open(filename, "w")
        setprint(fp, fp2, "# - - parameters - - - - - - - - - - - ")

        for (i, printlist_i) in enumerate(printlists)
            setprint(fp, fp2, printlists_header[i])
            for name in printlist_i
                headstring, key = get_header(params_set, name)
                if headstring !== nothing
                    if name == "meas_calls"
                        print_measurementinfo(fp, fp2, key)
                    else
                        string = get_stringfromkey(key)
                        paramstring = headstring * "[\"$name\"] = " * string
                        setprint(fp, fp2, paramstring)
                    end
                end
            end
            setprint(fp, fp2, "\t" )
        end
 
        setprint(fp, fp2, "# - - - - - - - - - - - - - - - - - - -")
        close(fp)
        close(fp2)
        println("""
        # Your parameters were written in $filename
        # If you want to do the simulation with same parameters, 
        # Just do 
        # julia run.jl $filename
        """)
        flush(stdout)
        return nothing
    end

    function print_parameters(params_set::ParamSet, p)
        print_parameters_list(params_set, p)
        return nothing
    end

    function print_parameters(filename, params_set::ParamSet)
        print_parameters_list(params_set, filename=filename)
		return nothing
    end


    function print_parameters(p)
        println("# - - parameters - - - - - - - - - - - ")
        pdict, pnames = make_parametersdict(p)

        for param in pdict
            if typeof(param[2]) == String
                println("$(param[1]) = \"$(param[2])\"")
            else
                println("$(param[1]) = $(param[2])")
            end
        end

        println("# - - - - - - - - - - - - - - - - - - -")
        print_parameters_file(p)
        println("""
        # Your parameters were written in parameters_used.jl
        # If you want to do the simulation with same parameters, 
        # Just do 
        # julia run.jl parameters_used.jl
        """)
        flush(stdout)
        return nothing
    end

    function parameterloading(physical, meta, param_meta, sim, mc, dirac, meas, system)
        param_set = ParamSet(physical, meta, param_meta, sim, mc, dirac, meas, system)
        p = Params(param_set)

        print_parameters(param_set, p)
        return p
    end

    function parameterloading(param_set::ParamSet)
        p = Params(param_set)

        print_parameters(param_set, p)
        return p
    end
    
end