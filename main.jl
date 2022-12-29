include("GaugeAction.jl")
include("Metropolis.jl")
include("TopCharge.jl")
include("ContCharge.jl")
include("Metadynamics.jl")
include("pretty_plots.jl")
include("error_est.jl")

using Plots
using Statistics
using StatsBase
using FFTW
using LaTeXStrings
using JLD

### Physics Parameters ###
const N_s = 32 
const N_t = 32
const d = 2
const β = N_s*N_t/80
### Metadynamics Parameters ###
const ε = 0.2 #Metropolis step-size
const Q_max = 10
const Q_thresh = 8
const dq = 0.001
const w = 1e-4
const k = 1000
### Simulation Parameters ###
const N_therm = 50000
const N_sweeps = 1000000
const n_skip = 10

boundary_count = 0
acc = 0

start = zeros(N_s,N_t,d);  #Startconfiguration

#bias_potential = zeros(grid_ind(Float64(Q_max))); #Bias-Potential Grid
q_vals = range(-Q_max,Q_max-dq,length(bias_potential));

Q_top,Q_cont,Q2_weighted,bias_potential = metropolis_meta(start,bias_potential,false); #Build Bias-Potential
#Q_top,Q_cont,Q2_weighted,bias_potential = metropolis_meta(start,bias_potential,true); äSimulate with static Bias-Potential (or dynammic if preferred)

χ_top,δχ_top = bootstrap(Q2_weighted,1000,true)
χ_top = round(χ_top,digits=3)

pretty_plots(q_vals,bias_potential,Q_cont,Q_top)
