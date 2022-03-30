# NumCME.jl
An extensible toolkit for direct numerical solution of the Chemical Master Equation based on the [Finite State Projection](https://www.atmos.colostate.edu/~munsky/Papers/JChemPhys_124_044104.pdf) and related algorithms.
## Features 
This package aims to offer dynamic, fast, and customizable methods for direct numerical integration of the Chemical Master Equation (CME) in Julia. Currently, it has:
- Transient solution of the CME for time-homogeneous reaction rates/propensities as well as time-varying reaction rates/propensities using Finite State Projection and related variants.
- Dynamic state spaces: The truncated state space is adapted on-the-fly to remove redundant states with low probabilities and add more states to ensure the truncation error is within user-specified tolerance.
- Customizability: Users can choose how the dynamic state space is managed (by specifying parameters for existing `SpaceAdapter` subtypes or write their own `SpaceAdapter`) and how the reduced ODEs are solved (by choosing one among the multitude options offered by [`DifferentialEquations.jl`](https://github.com/SciML/DifferentialEquations.jl)). Advanced users can write their own dynamic state space management policy by subtyping `SpaceAdapter`.
- Sensitivity analysis: Compute partial derivatives of the FSP solution with respect to model parameters. Users do not need to write analytic expressions of the propensity's partial derivatives. Instead, the package applies existing tools from [`ModelingToolkit.jl`](https://github.com/SciML/ModelingToolkit.jl) and [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl) to generate those derivatives automatically.
- Interface to [`Catalyst.jl`](https://catalyst.sciml.ai/dev/) package to accept reaction systems defined using Catalyst DSL.
## Installation
This package can be installed using Julia's package management. For the last stable version from the General Registry,
```julia
import Pkg; Pkg.add("NumCME")
```
For the latest commit on this repository, 
```julia
import Pkg; Pkg.add("https://github.com/voduchuy/NumCME.jl")
```
## Acknowledgement
Thanks to [Kaan Öcal](https://github.com/kaandocal) for the helpful conversation during [UQ-Bio 2021](https://q-bio.org/wp/qbss/) about Julia and Catalyst.