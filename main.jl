include("GaugeAction.jl")
include("Metropolis.jl")
include("TopCharge.jl")
include("ContCharge.jl")
include("Metadynamics.jl")
include("Bootstrap.jl")

using Plots
using Statistics
using StatsBase
using FFTW

const N_s = 24   #Definiere Theorie-Parameter (Anzahl an Zeit-cuts, und Potenzialkonstante)
const N_t = 24
const beta = 7.2
const epsilon = 0.2
const Q_max = 15
const dq = 0.01
const w = 1e-4

const N_therm = 1000
const N_sweeps = 100000   #Anzahl an Metropolis-Sweeps in der Simulation
boundary_count = 0

start = zeros(N_s,N_t,2).*(2*pi);  #Startkonfiguration
bias_potential = zeros(grid_ind(Q_max,dq,Q_max)); #Bias-Potential Grid

configs_meta = metropolis_meta(start,N_s,N_t,beta,epsilon,Q_max,dq,w,N_therm,N_sweeps,bias_potential);
#configs = metropolis(N_sweeps,N_therm,N_s,N_t,beta,epsilon)

Q_top_meta = topological_charge(configs_meta);
Q_cont_meta = cont_charge_full(configs_meta);

q_vals = range(-Q_max,Q_max,length(bias_potential));
#bias_plot = plot(q_vals,bias_potential,show=true)