using Random

N = 32
physical["N"] = (N,N)
physical["β"] = N^2 / 80
physical["Ntherm"] = 1_000
physical["Nsweeps"] = 10_000
physical["initial"] = "cold"

update["ϵ_metro"] = 0.3
update["metro_multi_hit"] = 1
update["metro_target_acc"] = 0.5

dirac["operator"] = "Wilson"
dirac["mass"] = 1.0
dirac["BC"] = [1, 1]

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
    Dict{Any,Any}("methodname" => "Chiral_condensate", "measure_every" => 100),
    Dict{Any,Any}("methodname" => "Dirac_determinant", "measure_every" => 100),
]
                       
system["veryverbose"] = false 
system["randomseeds"] = [Xoshiro()]

system["logdir"] = "./logs/N$(N)x$(N)_beta$(physical["β"])/"
system["logfile"] = "test"

system["measure_basedir"] = "./measurements/N$(N)x$(N)_beta$(physical["β"])/"
system["measure_dir"] = "test"

