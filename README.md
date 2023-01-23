# Schwinger-Metadynamics
Variations of Metadynamics in Simulations of the Schwinger Model (1+1D QED)

## Quick Start
1. Choose parameters of your choice to build a Metapotential first (use: template_build.jl to see structure) and do:
```
julia "build.jl" "build_parameters.jl"
```

2. Use built Metapotential to measure observables in second run (use: template_sym.jl to see structure) and do:
```
julia -t 2 "run.jl" "sim_parameters.jl"
```
**! Make sure to use 2 or more CPU-Threads when using Parallel Tempering!** 

**! Make sure to use the same Physics and Metapotential parameters for build- and measurement-run !**

Measurements are outputted as .txt files in the chosen directories; you can use them to make plots as you wish.
