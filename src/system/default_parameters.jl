defaultmeasures = Vector{Dict}(undef, 3)

for i in 1:length(defaultmeasures)
    defaultmeasures[i] = Dict()
end

defaultmeasures[1]["methodname"] = "Meta_charge"
defaultmeasures[1]["measure_every"] = 1
defaultmeasures[2]["methodname"] = "Topological_charge"
defaultmeasures[2]["measure_every"] = 1
defaultmeasures[3]["methodname"] = "Plaquette"
defaultmeasures[3]["measure_every"] = 1

physical = Dict()
meta = Dict()
param_meta = Dict()
opes = Dict()
update = Dict()
dirac = Dict()
meas = Dict()
system = Dict()

physical["starting_Q"] = nothing
physical["instanton_enabled"] = false
physical["swap_every"] = nothing

meta["meta_enabled"] = false
meta["tempering_enabled"] = false
meta["numinstances"] = 1
# meta["tempering_heatbath"] = false
meta["parametric"] = false
meta["symmetric"] = true
meta["is_static"] = [false]
meta["well_tempered"] = false
meta["take_snapshot_every"] = nothing
meta["no_zero_instance"] = false
meta["stride"] = 1

opes["is_static"] = [false]
opes["symmetric"] = true
opes["opes_enabled"] = false
opes["explore"] = nothing
opes["barrier"] = nothing
opes["biasfactor"] = nothing
opes["stride"] = nothing
opes["sigma0"] = nothing
opes["adaptive_σ_stride"] = nothing
opes["σ_min"] = nothing
opes["opes_epsilon"] = nothing
opes["cutoff"] = nothing
opes["d_thresh"] = nothing
opes["no_Z"] = nothing
opes["fixed_σ"] = nothing

param_meta["lower_bounds"] = nothing
param_meta["upper_bounds"] = nothing
param_meta["batchsize"] = nothing
param_meta["testfun"] = nothing
param_meta["minimizer"] = nothing

update["update_method"] = "Metro"
update["metro_ϵ"] = 0.2
update["metro_multi_hit"] = 1
update["metro_target_acc"] = 0.7
update["hmc_integrator"] = "Leapfrog"
update["hmc_steps"] = 10
update["hmc_Δτ"] = 0.7
update["hmc_ϕ"] = π/2

dirac["operator"] = nothing
dirac["mass"] = nothing
dirac["BC"] = nothing

meas["meas_calls"] = defaultmeasures

system["write_state_every"] = 10000
system["veryverbose"] = false
system["randomseeds"] = [Random.Xoshiro()]
