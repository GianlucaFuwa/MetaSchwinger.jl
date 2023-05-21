using Random
N = 32
physical["N"] = (N, N)
physical["β"] = N^2 / 80

meta["meta_enabled"] = true
    meta["parametric"] = true
    meta["symmetric"] = true
    meta["is_static"] = [false]
    meta["CVlims"] = (-7, 7)
    meta["bin_width"] = 1e-2
    meta["w"] = 1e-3
    meta["k"] = 1000
    meta["well_tempered"] = true
    meta["ΔT"] = 10.0
#
param_meta["potential_parameters"] = [-0.0, 20.0, 1.1]
param_meta["lower_bounds"] = [-0.5, 15.0, 1.0]
param_meta["upper_bounds"] = [-0.0, 25.0, 1.2]
param_meta["batchsize"] = 1000
param_meta["testfun"] = "GADLTtest"
param_meta["minimizer"] = "BFGS"
#
sim["Ntherm"] = 20_000
sim["Nsweeps"] = 100_000
sim["initial"] = "cold"

mc["ϵ_metro"] = 0.2
mc["multi_hit"] = 1
mc["metro_target_acc"] = 0.7

meas["meas_calls"] = Dict[Dict{Any,Any}("methodname" => "Meta_charge", "measure_every" => 1),
                          #Dict{Any,Any}("methodname" => "Topological_charge", "measure_every" => 10),
                          Dict{Any,Any}("methodname" => "Action", "measure_every" => 1),
                          #Dict{Any,Any}("methodname" => "Plaquette", "measure_every" => 10),
                          #Dict{Any,Any}("methodname" => "Wilson_loop_x16", "measure_every" => 10),
                          #Dict{Any,Any}("methodname" => "Polyakov_loop", "measure_every" => 10)
                          ]
                       
system["veryverbose"] = false
system["randomseeds"] = [Random.Xoshiro()]

system["logdir"] = "./logs/N$(physical["N"])/beta$(physical["β"])/parametric_MetaD"
system["logfile"] = "CVlims$(meta["CVlims"])_dcv$(meta["bin_width"])_w$(meta["w"])_build.txt"

system["measure_basedir"] = "./measurements/N$(physical["N"])/beta$(physical["β"])/parametric_MetaD"
system["measure_dir"] = "CVlims$(meta["CVlims"])_dcv$(meta["bin_width"])_w$(meta["w"])_build"

system["savebias_dir"] = "./metapotentials/N$(physical["N"])/beta$(physical["β"])/parametric"
system["biasfile"] = "CVlims$(meta["CVlims"])_dcv$(meta["bin_width"])_w$(meta["w"])_build"