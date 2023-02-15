using Random

physical["N"] = (32,32)
physical["β"] = 12.8

meta["meta_enabled"] = false

mc["ϵ_metro"] = 0.14

sim["Ntherm"] = 50_000
sim["Nsweeps"] = 5_000_000
sim["initial"] = "cold"

system["randomseeds"] = [Random.Xoshiro(1206)]

system["meas_calls"] = Dict[Dict{Any,Any}("methodname" => "Continuous_charge","measure_every" => 10),
						    Dict{Any,Any}("methodname" => "Topological_charge","measure_every" => 10),
                            Dict{Any,Any}("methodname" => "Action","measure_every" => 10),
                            Dict{Any,Any}("methodname" => "Plaquette","measure_every" => 10),]

                            
system["logdir"] = "./logs"
system["logfile"] = "N$(physical["N"])_beta$(physical["β"])"
system["measure_dir"] = "./measurements/N$(physical["N"])_beta$(physical["β"])"