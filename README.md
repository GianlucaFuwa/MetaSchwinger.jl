# Schwinger-Metadynamics
Variations of Metadynamics in Simulations of the Schwinger Model (1+1D QED)

## Features
- Basic quenched U(1) Lattice Gauge Theory simulations
- Metadynamics (doi: 10.1073/pnas.202427399) simulations incl. Well-Tempered Metadynamics (doi: 10.1103/PhysRevLett.100.020603)
- Parallel Tempering (Geyer, C. J. (1991). Markov chain Monte Carlo maximum likelihood. Computing Science and Statistics: Proc. 23rd Symp. Interface, 156–163.) incl. the possibility to add a different Metapotential to each of the instances

## Quick Start
1. Choose parameters of your choice to build a Metapotential first (use: template/parameters_build.jl to see structure) and do:
```
julia src/build.jl templates/build_parameters.jl
```

2. Use built Metapotential to measure observables in second run (use: template/parameters_meas.jl to see structure) and do:
```
julia -t 2 src/run.jl templates/sim_parameters.jl
```
**! Make sure to specify the number of CPU-threads to be equal to or more than the amount of instances when using Parallel Tempering!** 

Measurements are outputted as .txt files in the chosen directories; you can use them to make plots as you wish.
