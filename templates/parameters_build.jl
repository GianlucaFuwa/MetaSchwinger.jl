using Random

N = 32
vel = 1_000
ds = 0.01
w = 20.0/ds / vel

physical["N"] = (N, N)
physical["β"] = N^2 / 80
physical["Ntherm"] = 10_000
physical["Nsweeps"] = 1_000_000
physical["initial"] = "cold"

meta["meta_enabled"] = true
meta["tempering_enabled"] = false
meta["parametric"] = false
meta["symmetric"] = true
meta["stride"] = 1
meta["is_static"] = [false]
meta["CVlims"] = (-7, 7)
meta["bin_width"] = ds
meta["w"] = w
meta["k"] = 100
meta["well_tempered"] = false
system["write_state_every"] = 10_000
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
    Dict{Any,Any}("methodname" => "Action", "measure_every" => 1),
]

system["veryverbose"] = false
system["randomseeds"] = [Random.Xoshiro(1206)]

system["logdir"] = "./logs/N$(physical["N"])/beta$(physical["β"])/Master/chapter7/ft_metad"
system["logfile"] = "ds$(ds)_w$(w)_vel$(vel)"

system["measure_basedir"] = "./measurements/N$(physical["N"])/beta$(physical["β"])/Master/chapter7/ft_metad"
system["measure_dir"] = "ds$(ds)_w$(w)_vel$(vel)"

system["savebias_dir"] = "./metapotentials/N$(physical["N"])/beta$(physical["β"])/Master/chapter7/ft_metad/ds$(ds)_w$(w)_vel$(vel)"
system["biasfile"] = "snap"
