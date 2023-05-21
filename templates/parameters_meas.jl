using Random
N = 32
physical["N"] = (N,N)
physical["β"] = N^2 / 80

meta["meta_enabled"] = true
meta["is_static"] = [true, true, true, true]
meta["well_tempered"] = false
meta["symmetric"] = true
meta["CVlims"] = (-7,7)
meta["bin_width"] = 1e-2
meta["w"] = 0.0
meta["k"] = 1000
#meta["ΔT"] = 5.0

sim["Ntherm"] = 20_000
sim["Nsweeps"] = 1_000_000
sim["initial"] = "cold"

sim["tempering_enabled"] = true
sim["Ninstances"] = 5
sim["swap_every"] = 1

mc["ϵ_metro"] = 0.20
mc["multi_hit"] = 1
mc["metro_target_acc"] = 0.7

dirac["Dirac_operator"] = nothing
dirac["mass"] = 1.0
dirac["BC"] = [1, 1]

meas["meas_calls"] = Dict[
    Dict{Any,Any}("methodname" => "Meta_charge", "measure_every" => 10),
    Dict{Any,Any}("methodname" => "Topological_charge", "measure_every" => 10),
    Dict{Any,Any}("methodname" => "Action", "measure_every" => 10),
    Dict{Any,Any}("methodname" => "Plaquette", "measure_every" => 10),
    #Dict{Any,Any}("methodname" => "Polyakov_loop","measure_every" => 10),
    #=
    Dict{Any,Any}(
    "methodname" => "Dirac_eigenvalues",
    "measure_every" => 100,
    "nev" => 10, # Number of eigenvalues per measurement
    "which" => :SM, # :SM => smallest magnitude, :LM => largest magnitude
    "maxiter" => 3000, # Max. number of Arnoldi iterations
    ),
    Dict{Any,Any}("methodname" => "Dirac_determinant", "measure_every" => 100),
    =#
]
                       
system["veryverbose"] = false 
system["randomseeds"] = [Xoshiro(), Xoshiro(), Xoshiro(), Xoshiro(), Xoshiro()]

system["logdir"] = "./logs/N$(physical["N"])/beta$(physical["β"])/Seminar"
system["logfile"] = "CVlims$(meta["CVlims"])_dcv$(meta["bin_width"])_w$(meta["w"])_5instance"

system["measure_basedir"] = "./measurements/N$(physical["N"])/beta$(physical["β"])/Seminar"
system["measure_dir"] = "CVlims$(meta["CVlims"])_dcv$(meta["bin_width"])_w$(meta["w"])_5instance"

system["savebias_dir"] = "./metapotentials/N$(physical["N"])/beta$(physical["β"])/Seminar"
system["biasfile"] = "CVlims$(meta["CVlims"])_dcv$(meta["bin_width"])_5instance"
system["usebiases"] = [
    "./metapotentials/N$(physical["N"])/beta$(physical["β"])/CVlims(-7, 7)_dcv0.01_10.txt",
    "./metapotentials/N$(physical["N"])/beta$(physical["β"])/CVlims(-7, 7)_dcv0.01_15.txt",
    "./metapotentials/N$(physical["N"])/beta$(physical["β"])/CVlims(-7, 7)_dcv0.01_18.txt",
    "./metapotentials/N$(physical["N"])/beta$(physical["β"])/CVlims(-7, 7)_dcv0.01.txt",
]

