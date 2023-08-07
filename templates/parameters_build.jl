using Random

N = 24
physical["N"] = (N, N)
physical["β"] = N^2 / 80
physical["Ntherm"] = 10_000
physical["Nsweeps"] = 100_000
physical["initial"] = "cold"

meta["meta_enabled"] = true
meta["parametric"] = true
meta["symmetric"] = true
meta["is_static"] = [false]
meta["CVlims"] = (-7, 7)
meta["bin_width"] = 1e-2
meta["w"] = 1e-3
meta["k"] = 1000
meta["well_tempered"] = false
# meta["ΔT"] = 10.0

param_meta["potential_parameters"] = [-0.2, 3.0, 1.2]
param_meta["lower_bounds"] = [-0.5, 0.0, 1.0]
param_meta["upper_bounds"] = [-0.0, 10.0, 1.5]
param_meta["batchsize"] = 1000
param_meta["testfun"] = "KStest"
param_meta["minimizer"] = "SAMIN"


update["ϵ_metro"] = 0.2
update["metro_multi_hit"] = 1
update["metro_target_acc"] = 0.7

meas["meas_calls"] = Dict[Dict{Any,Any}("methodname" => "Meta_charge", "measure_every" => 1),
                          #Dict{Any,Any}("methodname" => "Topological_charge", "measure_every" => 10),
                          Dict{Any,Any}("methodname" => "Action", "measure_every" => 1),
                          #Dict{Any,Any}("methodname" => "Plaquette", "measure_every" => 10),
                          #Dict{Any,Any}("methodname" => "Wilson_loop_x16", "measure_every" => 10),
                          #Dict{Any,Any}("methodname" => "Polyakov_loop", "measure_every" => 10)
                          ]

system["veryverbose"] = false
system["randomseeds"] = [Random.Xoshiro()]

system["logdir"] = "./logs/N$(N)x$(N)_beta$(physical["β"])/parametric_MetaD"
system["logfile"] = "CVlims$(meta["CVlims"])_dcv$(meta["bin_width"])_w$(meta["w"])_build.txt"

system["measure_basedir"] = "./measurements/N$(N)x$(N)_beta$(physical["β"])/parametric_MetaD"
system["measure_dir"] = "CVlims$(meta["CVlims"])_dcv$(meta["bin_width"])_w$(meta["w"])_build"

system["savebias_dir"] = "./metapotentials/N$(N)x$(N)_beta$(physical["β"])/parametric"
system["biasfile"] = "CVlims$(meta["CVlims"])_dcv$(meta["bin_width"])_w$(meta["w"])_build"
