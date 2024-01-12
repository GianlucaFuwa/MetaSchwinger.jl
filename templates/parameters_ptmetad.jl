using Random

N = 32
var = "PTMetaD"

physical["N"] = (N,N)
physical["β"] = N^2 / 80
physical["Ntherm"] = 10_000
physical["Nsweeps"] = 100_000_000
physical["initial"] = "cold"

meta["meta_enabled"] = true
meta["tempering_enabled"] = true
meta["numinstances"] = 2
meta["swap_every"] = 1
meta["is_static"] = [true, true]
meta["well_tempered"] = false
meta["symmetric"] = true
meta["CVlims"] = (-7, 7)
meta["bin_width"] = 1e-2
meta["w"] = 0.0
meta["k"] = 1000

update["ϵ_metro"] = 0.20
update["metro_multi_hit"] = 1
update["metro_target_acc"] = 0.7

meas["meas_calls"] = Dict[
    Dict{Any,Any}("methodname" => "Meta_charge", "measure_every" => 10),
    Dict{Any,Any}("methodname" => "Topological_charge", "measure_every" => 10),
    Dict{Any,Any}("methodname" => "Plaquette", "measure_every" => 10),
    Dict{Any,Any}("methodname" => "Wilson_loop", "measure_every" => 10, "LX" => 8, "LT" => 8),
]

str = "Proceedings_edit/N$(N)x$(N)_beta$(N^2/80)"

system["veryverbose"] = false
system["randomseeds"] = [Xoshiro(), Xoshiro()]

system["logdir"] = "./logs/$(str)"
system["logfile"] = var

system["measure_basedir"] = "./measurements/$(str)"
system["measure_dir"] = var

system["savebias_dir"] = "./metapotentials/$(str)"
system["biasfile"] = var
system["usebiases"] = [
    "./metapotentials/$(str)/regular.txt",
]
