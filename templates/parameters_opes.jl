using Random

sigma = 0.1

N = 32
physical["N"] = (N,N)
physical["β"] = N^2 / 80
physical["Ntherm"] = 10_000
physical["Nsweeps"] = 200_000
physical["initial"] = "cold"

meta["meta_enabled"] = false
opes["opes_enabled"] = true
meta["CVlims"] = (-7, 7)
opes["barrier"] = 40
opes["biasfactor"] = Inf
opes["stride"] = 1
opes["sigma0"] = sigma
opes["σ_min"] = 1e-6
opes["fixed_σ"] = false
opes["d_thresh"] = 1
opes["no_Z"] = false
system["write_state_every"] = 10000

update["ϵ_metro"] = 0.20
update["metro_multi_hit"] = 1
update["metro_target_acc"] = 0.7

meas["meas_calls"] = Dict[
    Dict{Any,Any}("methodname" => "Meta_charge", "measure_every" => 10),
    Dict{Any,Any}("methodname" => "Topological_charge", "measure_every" => 10),
    Dict{Any,Any}("methodname" => "Action", "measure_every" => 10),
]

system["veryverbose"] = false
system["randomseeds"] = [Xoshiro(1206)]

system["logdir"] = "./logs/N$(physical["N"])/beta$(physical["β"])/Master/chapter7"
system["logfile"] = "sigma0_$(sigma)"

system["measure_basedir"] = "./measurements/N$(physical["N"])/beta$(physical["β"])/Master/chapter7"
system["measure_dir"] = "sigma0_$(sigma)"

system["savebias_dir"] = "./metapotentials/N$(physical["N"])/beta$(physical["β"])/Master/chapter7/sigma0_$(sigma)"
system["kernelsfp"] = "kernels"
system["statefp"] = "states"
