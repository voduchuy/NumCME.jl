# NumCME.jl
An extensible toolkit for direct numerical solution of the Chemical Master Equation based on the [Finite State Projection](https://www.atmos.colostate.edu/~munsky/Papers/JChemPhys_124_044104.pdf) and related algorithms.
## Features 
This package aims to offer dynamic, fast, and customizable methods for direct numerical integration of the Chemical Master Equation (CME) in Julia. Currently, it has:
- Transient solution of the CME for time-homogeneous reaction rates/propensities as well as time-varying reaction rates/propensities using Finite State Projection and related variants.
- Dynamic state spaces: The truncated state space is adapted on-the-fly to remove redundant states with low probabilities and add more states to ensure the truncation error is within user-specified tolerance.
- Customizable: Users can choose how the dynamic state space is managed (by specifying parameters for existing `SpaceAdapter` subtypes or write their own `SpaceAdapter`) and how the reduced ODEs are solved (by choosing one among the multitude options offered by [`DifferentialEquations.jl`](https://github.com/SciML/DifferentialEquations.jl)).
- Sensitivity analysis: Compute partial derivatives of the FSP solution with respect to model parameters. Users do not need to write analytic expressions of the propensity's partial derivatives. Instead, the package applies existing tools from [`ModelingToolkit.jl`](https://github.com/SciML/ModelingToolkit.jl) and [`ForwardDiff.jl`](https://github.com/JuliaDiff/ForwardDiff.jl) to generate those derivatives automatically.
- Advanced users can write their own dynamic state space management policy by subtyping `SpaceAdapter`.
## Installation
This package can be installed using Julia's package management. For the last stable version from the General Registry,
```julia
import Pkg; Pkg.add("NumCME")
```
For the latest commit on this repository, 
```julia
import Pkg; Pkg.add("https://github.com/voduchuy/NumCME.jl")
```
## TODO
- [ ] Tutorials.
- [x] User-friendly model input by integration with `Catalyst.jl`.
- [ ] Implement more state space adaptation policy:
    - [Sliding windows](https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-4-42).
    - [SSA-driven state space](https://doi.org/10.1016/j.mbs.2015.08.010).
    - [Multi-finite buffer](http://gila.bioe.uic.edu/lab/papers/2016/ACME-CaoTerebusLiang-2016.pdf).
- [ ] Implement [methods](https://doi.org/10.1063/1.4994917) based on Tensor-train approximations.
- [ ] [Stationary FSP](https://pubmed.ncbi.nlm.nih.gov/29055349/).
- [ ] User-friendly support for computing the [Fisher Information Matrix](https://doi.org/10.1371/journal.pcbi.1006365) and optimal experiment design.
## Other CME packages
The following packages are not written by Huy Vo but students and researchers interested in direct solution of the CME can additionally look at:
- [`FiniteStateProjection.jl`](https://github.com/kaandocal/FiniteStateProjection.jl) written by Kaan Öcal offers methods to convert a `Catalyst.jl` reaction network into a `ModelingToolkitize`'d FSP system (on a hyper-rectangular truncated state space) that can be solved using `DifferentialEquations.jl`. The package also automatically detects subnetworks with mass conservation (e.g. `A <->B` reactions) to reduce state dimension.
