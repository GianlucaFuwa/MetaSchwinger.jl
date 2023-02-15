using Random

physical["N"] = (104,104)
physical["β"] = 135.2

meta["Qmax"] = (-6,6)
meta["δq"] = 1e-4
meta["Qthr"] = (-6+meta["δq"],6-meta["δq"])
meta["w"] = 5e-4
meta["k"] = 1000

mc["ϵ_metro"] = 0.055

sim["Ntherm"] = 50_000
sim["Nsweeps"] = 1_000_000
sim["initial"] = "cold"
sim["parallel_tempering"] = true
sim["swap_every"] = 1

system["randomseeds"] = [Random.Xoshiro(1206),Random.Xoshiro(2805)]
                            
system["logdir"] = "./logs"
system["logfile"] = "N$(physical["N"])_beta$(physical["β"])_Qmax$(meta["Qmax"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_w$(meta["w"])_k$(meta["k"])_build.txt"
system["measure_dir"] = "./measurements/N$(physical["N"])_beta$(physical["β"])_Qmax$(meta["Qmax"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_w$(meta["w"])_k$(meta["k"])_build"
system["savebias_dir"] = "./metapotentials"
system["biasfile"] = "./N$(physical["N"])_beta$(physical["β"])_Qmax$(meta["Qmax"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_w$(meta["w"])_k$(meta["k"])_build.txt"
