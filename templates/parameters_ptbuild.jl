using Random
N = 32
physical["N"] = (N, N)
physical["β"] = N^2 / 80

meta["meta_enabled"] = true
    meta["parametric"] = false
    meta["symmetric"] = true
    meta["is_static"] = [false, false, false, false, false, false, false, false]
    meta["CVlims"] = (-7, 7)
    meta["bin_width"] = 1e-2
    meta["w"] = 2e-3
    meta["k"] = 1000
    meta["well_tempered"] = false
    #meta["ΔT"] = 5.0
    meta["take_snapshot_every"] = 10_000 * 8
#=
param_meta["potential_parameters"] = [-0.25, 19.5, 1.044]
param_meta["lower_bounds"] = [-0.26, 15.0, 1.0]
param_meta["upper_bounds"] = [-0.245, 25.0, 1.1]
param_meta["batchsize"] = 1000
param_meta["testfun"] = "GADLTtest"
param_meta["minimizer"] = "SAMIN"
=#
sim["Ntherm"] = 20_000
sim["Nsweeps"] = 250_000 * 8
sim["initial"] = "cold"
sim["starting_Q"] = [0, 1, 2, 3, 4, 5, 6, 7]
sim["tempering_enabled"] = true
    sim["Ninstances"] = 1
    sim["swap_every"] = nothing
    sim["no_zero_instance"] = true

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
system["randomseeds"] = [
    Random.Xoshiro(),
    Random.Xoshiro(),
    Random.Xoshiro(),
    Random.Xoshiro(),
    Random.Xoshiro(),
    Random.Xoshiro(),
    Random.Xoshiro(),
    Random.Xoshiro(),
]

system["logdir"] = "./logs/N$(physical["N"])/beta$(physical["β"])/Seminar"
system["logfile"] = "CVlims$(meta["CVlims"])_dcv$(meta["bin_width"])_w$(meta["w"])_normal_build"

system["measure_basedir"] = "./measurements/N$(physical["N"])/beta$(physical["β"])/Seminar"
system["measure_dir"] = "CVlims$(meta["CVlims"])_dcv$(meta["bin_width"])_w$(meta["w"])_normal_build"

system["savebias_dir"] = "./metapotentials/N$(physical["N"])/beta$(physical["β"])/Seminar"
system["biasfile"] = "CVlims$(meta["CVlims"])_dcv$(meta["bin_width"])_w$(meta["w"])_normal_build"
system["usebiases"] = [
    nothing,
]
system["snapshot_dir"] = "CVlims$(meta["CVlims"])_dcv$(meta["bin_width"])_w$(meta["w"])_normal_build_snapshots"