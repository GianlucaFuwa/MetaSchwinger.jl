using Random

N = 32
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
meta["bin_width"] = 0.01
meta["w"] = 0.02
meta["k"] = 1000
meta["well_tempered"] = false
system["write_state_every"] = 1000
# meta["ΔT"] = 10.0

param_meta["potential_parameters"] = [-0.2, 20.0, 1.1]
param_meta["lower_bounds"] = [-0.5, 15.0, 1.0]
param_meta["upper_bounds"] = [-0.0, 25.0, 1.5]
param_meta["batchsize"] = 1000
param_meta["testfun"] = "KStest"
param_meta["minimizer"] = "SAMIN"

update["ϵ_metro"] = 0.2
update["metro_multi_hit"] = 100
update["metro_target_acc"] = 0.5

meas["meas_calls"] = Dict[
    Dict{Any,Any}("methodname" => "Meta_charge", "measure_every" => 1),
    Dict{Any,Any}("methodname" => "Topological_charge", "measure_every" => 10),
    Dict{Any,Any}("methodname" => "Action", "measure_every" => 1),
]

system["veryverbose"] = false
system["randomseeds"] = [Random.Xoshiro()]

system["logdir"] = "./logs/N$(physical["N"])/beta$(physical["β"])/parametric_build1000_multi100"
system["logfile"] = "snap"

system["measure_basedir"] = "./measurements/N$(physical["N"])/beta$(physical["β"])"
system["measure_dir"] = "parametric_build1000_multi100"

system["savebias_dir"] = "./metapotentials/N$(physical["N"])/beta$(physical["β"])"
system["biasfile"] = "parametric_build1000_multi100"
