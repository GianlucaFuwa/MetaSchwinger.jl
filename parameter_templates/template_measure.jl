using Random

physical["N"] = (48,48)
physical["β"] = 28.8

meta["Qmax"] = (-6.5,6.5)
meta["Qthr"] = (-5.8,5.8)
meta["δq"] = 1e-3
meta["w"] = 1e-4
meta["k"] = 100

sim["ϵ"] = 0.2
sim["Ntherm"] = 25_000
sim["Nsweeps"] = 1_000_000
sim["insta_every"] = typemax(Int)
sim["meta_runtype"] = "dynamic"
sim["weightmode"] = "from_mean"
sim["initial"] = "cold"
sim["parallel_tempering"] = true
sim["swap_every"] = 100

system["randomseed"] = Random.Xoshiro(
    0x0c59f9d0d6a1e02b,
    0x9929c1fbac5d828f,
    0xf384432ed102e01f,
    0xd80d9f6f25af5c20)
system["meas_calls"] = Dict[Dict{Any,Any}("methodname" => "Continuous_charge","measure_every" => 10),
						    Dict{Any,Any}("methodname" => "Topological_charge","measure_every" => 10),
                            Dict{Any,Any}("methodname" => "Action","measure_every" => 10),
                            Dict{Any,Any}("methodname" => "Topological_susceptibility","measure_every" => 10)]

                            
system["logdir"] = "./logs"
system["logfile"] = "N$(physical["N"])_beta$(physical["β"])_Qmax$(meta["Qmax"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_w$(meta["w"])_k$(meta["k"])_temper$(sim["swap_every"]).txt"
system["measure_dir"] = "./measurements/N$(physical["N"])_beta$(physical["β"])_Qmax$(meta["Qmax"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_w$(meta["w"])_k$(meta["k"])_temper$(sim["swap_every"])"
system["savebias_dir"] = "./metapotentials"
system["biasfile"] = "N$(physical["N"])_beta$(physical["β"])_Qmax$(meta["Qmax"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_w$(meta["w"])_k$(meta["k"])_temper$(sim["swap_every"]).txt"
system["usebias"] = "./metapotentials/N$(physical["N"])_beta$(physical["β"])_Qmax$(meta["Qmax"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_w$(meta["w"])_k$(meta["k"])_build.txt"
