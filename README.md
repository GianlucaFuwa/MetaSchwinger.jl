# MetaSchwinger.jl
Variations of Metadynamics in Simulations of the Schwinger Model (1+1D QED)

## Features
- Basic quenched Schwinger Model simulations
- [Metadynamics](https://www.researchgate.net/publication/224908601_Metadynamics_A_method_to_simulate_rare_events_and_reconstruct_the_free_energy_in_biophysics_chemistry_and_material_science)
- [PT-MetaD](https://arxiv.org/abs/2307.04742)

## Quick Start
1. (Optional) Choose parameters of your choice to build a Metapotential first (use: template/parameters_build.jl to see structure) and do:
```
julia src/build.jl templates/build_parameters.jl
```

2. Use Metapotential to measure observables in second run (use: template/parameters_meas.jl to see structure) and do:
```
julia -t 2 src/run.jl templates/sim_parameters.jl
```
**! Make sure to specify the number of CPU-threads to be equal to or more than the amount of instances when using Parallel Tempering!** 

Measurements are outputted as .txt files in the chosen directories; you can use them to make plots as you wish.
