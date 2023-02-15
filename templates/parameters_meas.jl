using Random

physical["N"] = (104,104)
physical["β"] = 135.2

meta["Qmax"] = (-6,6)
meta["δq"] = 5e-4
meta["Qthr"] = (-6+meta["δq"],6-meta["δq"])
meta["w"] = 1e-5
meta["k"] = 1000

mc["ϵ_metro"] = 0.055

sim["Ntherm"] = 50_000
sim["Nsweeps"] = 5_000_000
sim["initial"] = "cold"
sim["parallel_tempering"] = true
sim["swap_every"] = 1

system["randomseeds"] = [Random.Xoshiro(1206),Random.Xoshiro(2805)]

system["meas_calls"] = Dict[Dict{Any,Any}("methodname" => "Continuous_charge","measure_every" => 10),
						    Dict{Any,Any}("methodname" => "Topological_charge","measure_every" => 10),
                            Dict{Any,Any}("methodname" => "Action","measure_every" => 10),
                            Dict{Any,Any}("methodname" => "Plaquette","measure_every" => 10),]

                            
system["logdir"] = "./logs"
system["logfile"] = "N$(physical["N"])_beta$(physical["β"])_Qmax$(meta["Qmax"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_w$(meta["w"])_k$(meta["k"])_TEMPER$(sim["swap_every"])_HIGHSTAT.txt"
system["measure_dir"] = "./measurements/N$(physical["N"])_beta$(physical["β"])_Qmax$(meta["Qmax"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_w$(meta["w"])_k$(meta["k"])_TEMPER$(sim["swap_every"])_HIGHSTAT"
system["savebias_dir"] = "./metapotentials"
system["biasfile"] = "N$(physical["N"])_beta$(physical["β"])_Qmax$(meta["Qmax"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_w$(meta["w"])_k$(meta["k"])_TEMPER$(sim["swap_every"])_HIGHSTAT.txt"
system["usebias"] = "./metapotentials/N$(physical["N"])_beta$(physical["β"])_Qmax$(meta["Qmax"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_w0.0001_k1000_build.txt"
