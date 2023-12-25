using Random

N = 40
physical["N"] = (N,N)
physical["β"] = N^2 / 80
physical["Ntherm"] = 10_000
physical["Nsweeps"] = 100_000
physical["initial"] = "cold"

opes["opes_enabled"] = true
meta["CVlims"] = (-6, 6)
opes["barrier"] = 40
opes["biasfactor"] = 40
opes["stride"] = 1
opes["sigma0"] = 0.0
opes["σ_min"] = 1e-6
opes["fixed_σ"] = false
opes["d_thresh"] = 1
opes["no_Z"] = false
system["write_state_every"] = 1000

update["ϵ_metro"] = 0.20
update["metro_multi_hit"] = 1
update["metro_target_acc"] = 0.7

meas["meas_calls"] = Dict[
    Dict{Any,Any}("methodname" => "Meta_charge", "measure_every" => 10),
    Dict{Any,Any}("methodname" => "Topological_charge", "measure_every" => 10),
    Dict{Any,Any}("methodname" => "Action", "measure_every" => 10),
    Dict{Any,Any}("methodname" => "Plaquette", "measure_every" => 10),
]

system["veryverbose"] = false
system["randomseeds"] = [Xoshiro(1206)]

system["logdir"] = "./logs/N$(physical["N"])/beta$(physical["β"])"
system["logfile"] = "test_opes"

system["measure_basedir"] = "./measurements/N$(physical["N"])/beta$(physical["β"])"
system["measure_dir"] = "test_opes"

system["savebias_dir"] = "./metapotentials/N$(physical["N"])/beta$(physical["β"])/test_opes_kernels"
system["kernelsfp"] = "snapshot"
system["statefp"] = "test_opes"
