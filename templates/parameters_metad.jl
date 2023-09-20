using Random

N = 32
physical["N"] = (N,N)
physical["β"] = N^2 / 80
physical["Ntherm"] = 10_000
physical["Nsweeps"] = 1_000_000
physical["initial"] = "cold"

meta["meta_enabled"] = true
meta["tempering_enabled"] = true
meta["numinstances"] = 3
meta["swap_every"] = 1
meta["is_static"] = [true, true, true]
meta["well_tempered"] = false
meta["symmetric"] = true
meta["CVlims"] = (-7,7)
meta["bin_width"] = 1e-2
meta["w"] = 0.0
meta["k"] = 1000

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
system["randomseeds"] = [Xoshiro(1206), Xoshiro(2805), Xoshiro(0505)]

system["logdir"] = "./logs/N$(physical["N"])/beta$(physical["β"])/Seminar"
system["logfile"] = "lukas_0"

system["measure_basedir"] = "./measurements/N$(physical["N"])/beta$(physical["β"])/Seminar"
system["measure_dir"] = "lukas_0"

system["savebias_dir"] = "./metapotentials/N$(physical["N"])/beta$(physical["β"])/Seminar"
system["biasfile"] = "lukas_0"
system["usebiases"] = [
    "./metapotentials/N$(physical["N"])/beta$(physical["β"])/CVlims(-7, 7)_dcv0.01_notflat.txt",
    "./metapotentials/N$(physical["N"])/beta$(physical["β"])/CVlims(-7, 7)_dcv0.01_notflat.txt",
]
