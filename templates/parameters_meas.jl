using Random
N = 28
physical["N"] = (N,N)
physical["β"] = N^2 / 80

meta["meta_enabled"] = true
    meta["is_static"] = [false]
    meta["well_tempered"] = false
    meta["symmetric"] = true
    meta["CVlims"] = (-7,7)
    meta["bin_width"] = 1e-2
    meta["w"] = 1e-4
    meta["k"] = 100000
    #meta["ΔT"] = 5.0


sim["Ntherm"] = 20_000
sim["Nsweeps"] = 100_000_000
sim["initial"] = "cold"

sim["tempering_enabled"] = true
    sim["Ninstances"] = 2
    sim["swap_every"] = 1

mc["ϵ_metro"] = 0.20

meas["meas_calls"] = Dict[Dict{Any,Any}("methodname" => "Meta_charge","measure_every" => 10),
                          Dict{Any,Any}("methodname" => "Topological_charge","measure_every" => 10),
                          Dict{Any,Any}("methodname" => "Action","measure_every" => 10),
                          Dict{Any,Any}("methodname" => "Plaquette","measure_every" => 10),
                          #Dict{Any,Any}("methodname" => "Polyakov_loop","measure_every" => 10),
                          ]
                       
system["veryverbose"] = false 
system["randomseeds"] = [Random.Xoshiro(1206), Random.Xoshiro()]

system["logdir"] = "./logs/N$(physical["N"])/beta$(physical["β"])/Proceedings"
system["logfile"] = "CVlims$(meta["CVlims"])_dcv$(meta["bin_width"])_w$(meta["w"])_meas_PT.txt"

system["measure_basedir"] = "./measurements/N$(physical["N"])/beta$(physical["β"])/Proceedings"
system["measure_dir"] = "CVlims$(meta["CVlims"])_dcv$(meta["bin_width"])_w$(meta["w"])_meas_PT"

system["savebias_dir"] = "./metapotentials/N$(physical["N"])/beta$(physical["β"])"
system["biasfile"] = "CVlims$(meta["CVlims"])_dcv$(meta["bin_width"])_meas_PT"
system["usebiases"] = ["./metapotentials/N$(physical["N"])/beta$(physical["β"])/CVlims(-7, 7)_dcv0.01.txt"]

