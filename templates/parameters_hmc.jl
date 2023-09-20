using Random

N = 32
physical["N"] = (N,N)
physical["β"] = N^2 / 80
physical["Ntherm"] = 1_00
physical["Nsweeps"] = 100
physical["initial"] = "cold"

meta["meta_enabled"] = true
meta["tempering_enabled"] = false
meta["numinstances"] = 1
meta["swap_every"] = 1
meta["is_static"] = [true]
meta["well_tempered"] = false
meta["symmetric"] = true
meta["CVlims"] = (-7,7)
meta["bin_width"] = 1e-2
meta["w"] = 0.0
meta["k"] = 1000

update["update_method"] = "hmc"
update["hmc_integrator"] = "OMF4"
update["hmc_Δτ"] = 0.2
update["hmc_steps"] = 5
update["hmc_ϕ"] = π/2

meas["meas_calls"] = Dict[
    Dict{Any,Any}("methodname" => "Meta_charge", "measure_every" => 10),
    Dict{Any,Any}("methodname" => "Topological_charge", "measure_every" => 10),
    Dict{Any,Any}("methodname" => "Action", "measure_every" => 10),
    Dict{Any,Any}("methodname" => "Plaquette", "measure_every" => 10),
    Dict{Any,Any}("methodname" => "Polyakov_loop","measure_every" => 10),
    # Dict{Any,Any}(
    # "methodname" => "Dirac_eigenvalues",
    # "measure_every" => 100,
    # "nev" => 10, # Number of eigenvalues per measurement
    # "which" => :SM, # :SM => smallest magnitude, :LM => largest magnitude
    # "maxiter" => 3000, # Max. number of Arnoldi iterations
    # ),
    # Dict{Any,Any}("methodname" => "Chiral_condensate", "measure_every" => 100),
    # Dict{Any,Any}("methodname" => "Dirac_determinant", "measure_every" => 100),
]

system["veryverbose"] = false
system["randomseeds"] = [Xoshiro(1206)]

system["logdir"] = "./logs/N$(N)x$(N)_beta$(physical["β"])/Seminar"
system["logfile"] = "hmc_tau1_nofriction"

system["measure_basedir"] = "./measurements/N$(N)x$(N)_beta$(physical["β"])/Seminar"
system["measure_dir"] = "hmc_tau1_nofriction"

system["savebias_dir"] = "./metapotentials/N$(physical["N"])/beta$(physical["β"])/Seminar"
system["biasfile"] = "hmc_tau1_nofriction"
system["usebiases"] = [
    "./metapotentials/N$(physical["N"])/beta$(physical["β"])/CVlims(-7, 7)_dcv0.01_notflat.txt",
]
