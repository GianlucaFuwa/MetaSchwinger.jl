module System_parameters
	using Random
    using DelimitedFiles
	export Params
	
	printlist_physical = ["N","β"]

	printlist_meta = ["Qmax","Qthr","δq","w","k"]

	printlist_sim = ["ϵ","Ntherm","Nsweeps","insta_every","meta_runtype","weightmode","initial","parallel_tempering","swap_every"]

    printlist_system = ["randomseed","meas_calls","logdir","logfile","loadfile",
    "measure_dir","savebias_dir","biasfile","usebias","weightfile"]

    const printlists = [printlist_physical,printlist_meta,printlist_sim,printlist_system]
    const printlists_header = ["# Physical Settings ","# Metadynamics Settings",
    "# Simulation Settings","# System Settings"]

	defaultmeasures = Array{Dict,1}(undef,2)
	for i=1:length(defaultmeasures)
		defaultmeasures[i] = Dict()
	end
	defaultmeasures[1]["methodname"] = "Continuous_charge"
    defaultmeasures[1]["measure_every"] = 1
    defaultmeasures[2]["methodname"] = "Topological_charge"
    defaultmeasures[2]["measure_every"] = 1

	physical = Dict()
	meta = Dict()
	sim = Dict()
	system = Dict()

    sim["insta_every"] = typemax(Int)
    sim["meta_runtype"] = "dynamic"
    sim["weightmode"] = "from_mean"
    sim["parallel_tempering"] = false
    sim["swap_every"] = 10

    system["randomseed"] = Random.Xoshiro(
        0x0c59f9d0d6a1e02b,
        0x9929c1fbac5d828f,
        0xf384432ed102e01f,
        0xd80d9f6f25af5c20)
	system["meas_calls"] = defaultmeasures
	
	mutable struct Params_set
	physical::Dict
	meta::Dict
	sim::Dict
	system::Dict
	
        function Params_set(physical,meta,sim,system)
            return new(physical,meta,sim,system)
        end
    end
	
    function make_parameters(physical,meta,sim,system)
        return Params_set(physical,meta,sim,system)
    end
    
	struct Params
		N::Tuple
		β::Float64

		Qmax::NTuple{2,Float64}
		Qthr::NTuple{2,Float64}
		δq::Float64
		w::Float64
		k::Float64

		ϵ::Float64
		Ntherm::Int64
		Nsweeps::Int64
		insta_every::Int64
		meta_runtype::String
        weightmode::String
		initial::String
        parallel_tempering::Bool
        swap_every::Int64

		randomseed::Xoshiro
		meas_calls::Array{Dict,1}
		logdir::String
		logfile::String
		loadfile::IOStream
		measure_dir::String
        measure_dir_secondary::String
		savebias_dir::String
		biasfile::String
		usebias::Array{Float64,1}
        weightfile::String

		function Params(physical,meta,sim,system)
			N = physical["N"]
			β = physical["β"]

			Qmax = meta["Qmax"]
			Qthr = meta["Qthr"]
			δq = meta["δq"]
			w = meta["w"]
			k = meta["k"]
			
			ϵ = sim["ϵ"]
			Ntherm = sim["Ntherm"]
			Nsweeps = sim["Nsweeps"]
			insta_every = sim["insta_every"]
            if sim["meta_runtype"] == "static" || sim["meta_runtype"] == "dynamic"
			    meta_runtype = sim["meta_runtype"]
            else 
                error("meta_runtype can only be either \"static\" or \"dynamic\".")
            end
            if sim["weightmode"] == "from_mean" || sim["weightmode"] == "from_current"
                weightmode = sim["weightmode"]
            else 
                error("weightmode can only be either \"from_mean\" or \"from_current\".")
            end
			initial = sim["initial"]
            parallel_tempering = sim["parallel_tempering"]
            swap_every = sim["swap_every"]

			randomseed = system["randomseed"]
			meas_calls = system["meas_calls"]
			logdir = system["logdir"]
			if isdir(logdir) == false
				mkdir(logdir)
			end
			logfile = pwd()*"/"*logdir*"/"*system["logfile"]
			loadfile = open(logfile,"a")
			measure_dir = system["measure_dir"]
            if parallel_tempering
                measure_dir_secondary = system["measure_dir"]*"_secondary"
                weightfile = pwd()*"/"*measure_dir_secondary*"/Weights.txt"
            else 
                measure_dir_secondary = ""
                weightfile = pwd()*"/"*measure_dir*"/Weights.txt"
            end
			if isdir(measure_dir) == false
				mkpath(measure_dir)
            end
            if isdir(measure_dir_secondary) == false && measure_dir_secondary !== ""
                mkpath(measure_dir_secondary)
			end
			savebias_dir = system["savebias_dir"]
			if isdir(savebias_dir) == false
				mkpath(savebias_dir)
			end
			biasfile = pwd()*"/"*savebias_dir*"/"*system["biasfile"]
			if haskey(system,"usebias")
                file = system["usebias"]
				usebias = readdlm(file,Float64)
                @assert length(usebias[:,2]) == round(Int,((Qmax[2]-Qmax[1])/δq),RoundNearestTiesAway)+1 "Length of given Metapotential doesn't match Meta-parameters"
                usebias = usebias[:,2]
			else 
				usebias = zeros(round(Int,(Qmax[2]-Qmax[1])/δq,RoundNearestTiesAway)+1)
			end
			return new(
				N,β,Qmax,Qthr,δq,w,k,ϵ,Ntherm,Nsweeps,insta_every,meta_runtype,
                weightmode,initial,parallel_tempering,swap_every,randomseed,meas_calls,logdir,
				logfile,loadfile,measure_dir,measure_dir_secondary,savebias_dir,biasfile,usebias,weightfile)
		end

		function Params(params_set::Params_set)
			return Params(params_set.physical,params_set.meta,params_set.sim,params_set.system)
		end
	end
    
    function set_params(dict,string)
        if haskey(dict,string)
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
		return pdict,pnames
	end

	function print_parameters_file(p)
		filename = p.logfile*"_parameters.jl"
		fp = open(filename,"w")
		println(fp,"# - - parameters - - - - - - - - ")
		pdict, = make_parametersdict(p)
		for param in pdict
			if typeof(param[2]) == String
				println(fp,"$(param[1]) = \"$(param[2])\"")
			else
				println(fp,"$(param[1]) = $(param[2])")
			end
		end
		println(fp,"# - - - - - - - - - - - - - - - -")
		close(fp)
	end
	
    function setprint(fp,fp2,string)
        println(fp,string)
        if fp2 !== nothing
            println(fp2,string)
        end
        println(string)
    end
    
    function get_stringfromkey(key)
        if typeof(key) == String
            string = "\"$key\""
        else
            string = "$key"
		end
	end

	function get_header(params_set::Params_set,inputname)
		if haskey(params_set.system,inputname)
            return "system",params_set.system[inputname]
        elseif haskey(params_set.physical,inputname)
            return "physical",params_set.physical[inputname]
        elseif haskey(params_set.meta,inputname)
            return "meta",params_set.meta[inputname]
        elseif haskey(params_set.sim,inputname)
            return "sim",params_set.sim[inputname]
        else
            return nothing,nothing
        end
    end

	function print_measurementinfo(fp,fp2,key)
        string = "measurement[\"measurement_methods\"] = Dict[ "
        setprint(fp,fp2,string)
        for (i,data) in enumerate(key)
            string = "  Dict{Any,Any}(\"methodname\" => \"$(data["methodname"])\","
            setprint(fp,fp2,string)
            count = 0
            for (name,key_i) in data
                if name != "methodname"
                    count += 1
                    paramstring = get_stringfromkey(key_i)
                    string = "    \"$(name)\" => "*paramstring
                    if count != length(data) -1
                        string *= ","
                    end
                    if name == "fermiontype" && key_i === nothing
                    else
                        setprint(fp,fp2,string)
                    end
                    
                end
            end
            string = "  )"
            if i != length(key)
                string *= ","
            end 
            setprint(fp,fp2,string)
            
        end
        string = "]"
        setprint(fp,fp2,string)
    end

    function print_parameters_list(params_set::Params_set,p=nothing;filename=nothing)
        if filename === nothing 
            @assert p !== nothing "wrong input!"

            filename = p.logfile*"_parameters.jl"
            fp2 = p.loadfile
        else
            fp2 = nothing
        end
        fp = open(filename,"w")
        
        setprint(fp,fp2,"# - - parameters - - - - - - - - - - - ")
        for (i,printlist_i) in enumerate(printlists)
            setprint(fp,fp2,printlists_header[i] )
            for name in printlist_i
                headstring,key = get_header(params_set,name)
                if headstring !== nothing
                    if name == "measurement_methods"
                        print_measurementinfo(fp,fp2,key)
                    else
                        string = get_stringfromkey(key)
                        paramstring = headstring*"[\"$name\"] = "*string
                        setprint(fp,fp2,paramstring)
                    end
                end
            end
            setprint(fp,fp2,"\t" )
        end
 
        setprint(fp,fp2,"# - - - - - - - - - - - - - - - - - - -")

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

    function print_parameters(params_set::Params_set,p)
        print_parameters_list(params_set,p)
        return nothing
    end

    function print_parameters(filename,params_set::Params_set)
        print_parameters_list(params_set,filename=filename)
		return nothing
    end


    function print_parameters(p)
        println("# - - parameters - - - - - - - - - - - ")
        
        pdict,pnames = make_parametersdict(p)
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

    function parameterloading(physical,meta,sim,system)
        param_set = Params_set(physical,meta,sim,system)
        p = Params(param_set)

        print_parameters(param_set,p)
        return p
    end

    function parameterloading(param_set::Params_set)
        p = Params(param_set)

        print_parameters(param_set,p)
        return p
    end
    
end
