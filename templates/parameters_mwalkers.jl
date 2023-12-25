using Random
N = 32
physical["N"] = (N, N)
physical["β"] = N^2 / 80
physical["Ntherm"] = 10_000
physical["Nsweeps"] = 100_000
physical["initial"] = "cold"
physical["starting_Q"] = nothing

# meta["meta_enabled"] = true
# meta["parametric"] = false
meta["tempering_enabled"] = false
meta["numinstances"] = 1
# meta["swap_every"] = nothing
# meta["no_zero_instance"] = true
# meta["symmetric"] = true
# meta["is_static"] = [false for _ in 1:15]
# meta["bin_width"] = 0.01
# meta["w"] = 0.002
# meta["k"] = 1000
# meta["CVlims"] = (-7, 7)

meta["meta_enabled"] = false
opes["symmetric"] = true
opes["opes_enabled"] = true
opes["CVlims"] = (-7, 7)
opes["barrier"] = 40
opes["biasfactor"] = Inf
opes["stride"] = 1
opes["sigma0"] = 0.1
opes["σ_min"] = 1e-6
opes["fixed_σ"] = false
opes["d_thresh"] = 1
opes["no_Z"] = false
opes["stride"] = 1

system["write_state_every"] = 1000

# param_meta["potential_parameters"] = [-0.2, 20.0, 1.1]
# param_meta["lower_bounds"] = [-0.5, 15.0, 1.0]
# param_meta["upper_bounds"] = [-0.0, 25.0, 1.5]
# param_meta["batchsize"] = 100
# param_meta["testfun"] = "GADLTtest"
# param_meta["minimizer"] = "SAMIN"

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
system["randomseeds"] = [Xoshiro() for _ in 1:15]

system["logdir"] = "./logs/N$(physical["N"])/beta$(physical["β"])"
system["logfile"] = "test_opes1340_3"

system["measure_basedir"] = "./measurements/N$(physical["N"])/beta$(physical["β"])"
system["measure_dir"] = "test_opes1340_3"

system["savebias_dir"] = "./metapotentials/N$(physical["N"])/beta$(physical["β"])/test_opes1340_3"
system["biasfile"] = "snapshot"
system["kernelsfp"] = "kernels"
system["statefp"] = "states"
