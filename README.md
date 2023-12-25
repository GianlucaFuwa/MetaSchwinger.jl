# MetaSchwinger.jl
Variations of Metadynamics in Simulations of the Schwinger Model (1+1D QED)

## Features
- Basic quenched Schwinger Model simulations
- [Metadynamics](https://www.researchgate.net/publication/224908601_Metadynamics_A_method_to_simulate_rare_events_and_reconstruct_the_free_energy_in_biophysics_chemistry_and_material_science)
- [PT-MetaD](https://arxiv.org/abs/2307.04742)
- [OPES](https://pubs.acs.org/doi/abs/10.1021/acs.jpclett.0c00497)

## Quick Start
1. (Optional) Choose parameters of your choice to build a Metapotential first (use: template/parameters_build.jl to see structure) and do:
```
julia src/build.jl templates/build_parameters.jl
```

2. Use Metapotential to measure observables in second run (use: template/parameters_meas.jl to see structure) and do:
```
julia src/run.jl templates/sim_parameters.jl
```

Measurements are outputted as .txt files in the chosen directories; you can use them to make plots as you wish.
