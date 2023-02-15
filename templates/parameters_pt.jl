using Random

physical["N"] = (32,32)
physical["β"] = 12.8

meta["meta_enabled"] = true
meta["Qmax"] = (-6,6)
meta["δq"] = 5e-4
meta["Qthr"] = (-6+meta["δq"],6-meta["δq"])
meta["w"] = 1e-5
meta["k"] = 1000

mc["ϵ_metro"] = 0.14

sim["Ntherm"] = 1_000
sim["Nsweeps"] = 1_000
sim["initial"] = "cold"
sim["tempering_enabled"] = true
sim["Ninstances"] = 4
sim["swap_every"] = 1


meas["meas_calls"] = Dict[Dict{Any,Any}("methodname" => "Continuous_charge","measure_every" => 10),
                          Dict{Any,Any}("methodname" => "Topological_charge","measure_every" => 10),
                          Dict{Any,Any}("methodname" => "Action","measure_every" => 10),
                          Dict{Any,Any}("methodname" => "Plaquette","measure_every" => 10),]

system["veryverbose"] = false # NOT recommended for high statistic
system["randomseeds"] = [Random.Xoshiro(),Random.Xoshiro(),Random.Xoshiro(),Random.Xoshiro()]
system["logdir"] = "./logs/N$(physical["N"])/beta$(physical["β"])/PT1_4instances"
system["logfile"] = "Qmax$(meta["Qmax"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_w$(meta["w"])"
system["measure_dir"] = "./measurements/N$(physical["N"])/beta$(physical["β"])/Qmax$(meta["Qmax"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_w$(meta["w"])"
system["savebias_dir"] = "./metapotentials/N$(physical["N"])/beta$(physical["β"])/PT1_4instances"
system["biasfile"] = "Qmax$(meta["Qmax"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_w$(meta["w"])"
