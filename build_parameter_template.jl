using Random

physical["N"] = (48,48)
physical["β"] = 28.8

meta["Qmax"] = (-6.5,6.5)
meta["Qthr"] = (-5.8,5.8)
meta["δq"] = 1e-3
meta["w"] = 1e-4
meta["k"] = 100

sim["ϵ"] = 0.15
sim["Ntherm"] = 25_000
sim["Nsweeps"] = 1_000_000
sim["initial"] = "cold"
                            
system["logdir"] = "./logs"
system["logfile"] = "N$(physical["N"])_beta$(physical["β"])_Qmax$(meta["Qmax"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_w$(meta["w"])_k$(meta["k"])_build.txt"
system["measure_dir"] = "./measurements/N$(physical["N"])_beta$(physical["β"])_Qmax$(meta["Qmax"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_w$(meta["w"])_k$(meta["k"])_build"
system["savebias_dir"] = "./metapotentials"
system["biasfile"] = "N$(physical["N"])_beta$(physical["β"])_Qmax$(meta["Qmax"])_Qthr$(meta["Qthr"])_dq$(meta["δq"])_w$(meta["w"])_k$(meta["k"])_build.txt"
