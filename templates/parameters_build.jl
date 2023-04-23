using Random

physical["N"] = (8, 8)
physical["β"] = 0.8

meta["meta_enabled"] = true
meta["parametric"] = false
meta["symmetric"] = true
meta["is_static"] = [false]
meta["CVlims"] = (-7, 7)
meta["bin_width"] = 1e-2
meta["w"] = 1e-3
meta["k"] = 1000
meta["well_tempered"] = false
#meta["ΔT"] = 5.0
#=
param_meta["potential_parameters"] = [-0.0001, 15.0, 1.1]
param_meta["lower_bounds"] = [-0.5, 10.0, 1.0]
param_meta["upper_bounds"] = [0.0, 30.0, 1.2]
param_meta["batchsize"] = 1000
param_meta["testfun"] = "GADLTtest"
param_meta["minimizer"] = "SAMIN"
=#
sim["Ntherm"] = 20_000
sim["Nsweeps"] = 2_000_000
sim["initial"] = "cold"
sim["parallel_tempering"] = nothing
sim["swap_every"] = nothing

mc["ϵ_metro"] = 1.25
meas["meas_calls"] = Dict[Dict{Any,Any}("methodname" => "Meta_charge", "measure_every" => 10),
                          Dict{Any,Any}("methodname" => "Topological_charge", "measure_every" => 10),
                          Dict{Any,Any}("methodname" => "Action", "measure_every" => 10)]
                          #Dict{Any,Any}("methodname" => "Plaquette", "measure_every" => 10),
                          #Dict{Any,Any}("methodname" => "Wilson_loop_x16", "measure_every" => 10),
                          #Dict{Any,Any}("methodname" => "Polyakov_loop", "measure_every" => 10)]
                       
system["veryverbose"] = false
system["randomseeds"] = [Random.Xoshiro(1206)]

system["logdir"] = "./logs/N$(physical["N"])/beta$(physical["β"])/MetaD"
system["logfile"] = "CVlims$(meta["CVlims"])_dcv$(meta["bin_width"])_w$(meta["w"])_build.txt"

system["measure_basedir"] = "./measurements/N$(physical["N"])/beta$(physical["β"])/MetaD/"
system["measure_dir"] = "CVlims$(meta["CVlims"])_dcv$(meta["bin_width"])_w$(meta["w"])_build"

system["savebias_dir"] = "./metapotentials/N$(physical["N"])/beta$(physical["β"])"
system["biasfile"] = "CVlims$(meta["CVlims"])_dcv$(meta["bin_width"])_w$(meta["w"])_build"