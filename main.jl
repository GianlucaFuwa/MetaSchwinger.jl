include("GaugeAction.jl")
include("Metropolis.jl")
include("TopCharge.jl")
include("ContCharge.jl")
include("Metadynamics.jl")
include("pretty_plots.jl")

using Plots
using Statistics
using StatsBase
using FFTW
using LaTeXStrings

const N_s = 32   #Definiere Theorie-Parameter (Anzahl an Zeit-cuts, und Potenzialkonstante)
const N_t = 32
const d = 2
const β = N_s*N_t/80
const ε = 0.2
const Q_max = 15
const Q_thresh = 12
const dq = 0.001
const w = 1e-3
const k = 1000
const n_skip = 10

const N_therm = 50000
const N_sweeps = 2000000   #Anzahl an Metropolis-Sweeps in der Simulation
boundary_count = 0
acc = 0

start = zeros(N_s,N_t,d).*(2*pi);  #Startkonfiguration
bias_potential = zeros(grid_ind(Float64(Q_max))); #Bias-Potential Grid
q_vals = range(-Q_max,Q_max,length(bias_potential));

Q_top,Q_cont,Q2_weighted = metropolis_meta(start,bias_potential);
#configs = metropolis(N_sweeps,N_therm,N_s,N_t,beta,epsilon)

χ_top = sum(Q2_weighted)

pretty_plots(q_vals,bias_potential,Q_cont,Q_top)