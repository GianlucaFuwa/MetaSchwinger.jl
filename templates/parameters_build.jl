using Random

N = 32
physical["N"] = (N, N)
physical["β"] = N^2 / 80
physical["Ntherm"] = 10_000
physical["Nsweeps"] = 2_000_000
physical["initial"] = "cold"

meta["meta_enabled"] = true
meta["tempering_enabled"] = false
meta["parametric"] = false
meta["symmetric"] = true
meta["is_static"] = [false]
meta["CVlims"] = (-7, 7)
meta["bin_width"] = 0.01
meta["w"] = 0.002
meta["k"] = 1000
meta["well_tempered"] = false
system["write_state_every"] = 30000
# meta["ΔT"] = 10.0
opes["opes_enabled"] = false

# param_meta["potential_parameters"] = [-0.2, 3.0, 1.2]
# param_meta["lower_bounds"] = [-0.5, 0.0, 1.0]
# param_meta["upper_bounds"] = [-0.0, 10.0, 1.5]
# param_meta["batchsize"] = 1000
# param_meta["testfun"] = "KStest"
# param_meta["minimizer"] = "SAMIN"

update["ϵ_metro"] = 0.2
update["metro_multi_hit"] = 1
update["metro_target_acc"] = 0.7

meas["meas_calls"] = Dict[
    Dict{Any,Any}("methodname" => "Meta_charge", "measure_every" => 1),
    Dict{Any,Any}("methodname" => "Topological_charge", "measure_every" => 10),
    Dict{Any,Any}("methodname" => "Action", "measure_every" => 1),
]

system["veryverbose"] = false
system["randomseeds"] = [Random.Xoshiro()]

system["logdir"] = "./logs/N$(physical["N"])/beta$(physical["β"])"
system["logfile"] = "test_metad1"

system["measure_basedir"] = "./measurements/N$(physical["N"])/beta$(physical["β"])"
system["measure_dir"] = "test_metad1"

system["savebias_dir"] = "./metapotentials/N$(physical["N"])/beta$(physical["β"])/test_metad1"
system["biasfile"] = "snapshot"
