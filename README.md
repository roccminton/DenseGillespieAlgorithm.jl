# DenseGillespieAlgorithm

[![Build Status](https://github.com/roccminton/DenseGillespieAlgorithm.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/roccminton/DenseGillespieAlgorithm.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://roccminton.github.io/DenseGillespieAlgorithm.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://roccminton.github.io/DenseGillespieAlgorithm.jl/dev)

# DenseGillespieAlgorithm.jl

DenseGillespieAlgorithm.jl is a Julia package designed for performing exact stochastic simulations of systems with dense reaction networks. This package implements a modified Gillespie algorithm that is highly optimized for problems where the number of different species is high, but the number of reaction events stays small. It is ideal for simulating stochastic population dynamics, and other systems that require exact trajectory tracking in systems with a large number of interacting species.

## Features

- **Efficient Exact Stochastic Simulations**: Ideal for densely connected reaction networks.
- **Modified Gillespie Algorithm**: Optimized for dense problems, reducing the computational complexity.
- **Flexible**: Works for a variety of problem domains including chemical kinetics and population dynamics.

## Installation

### Install from GitHub

You can install the package directly from this GitHub repository:

```julia
using Pkg
Pkg.add(url="https://github.com/username/DenseGillespieAlgorithm.jl")
```
### Install from Julia

Once the package is registered in the official Julia package registry, you can install it via:

```julia
using Pkg
Pkg.add("DenseGillespieAlgorithm")
```

This will install the latest stable version of the package and all required dependencies.

## Basic Usage

The following minimal working example demonstrates the use of  **DenseGillespieAlgorithm.jl** for the simulation of a simple SIR model. It should be noted that the framework was not designed for such simple systems. However, the [JumpProcesses.jl](https://docs.sciml.ai/JumpProcesses/stable/) within the SciML ecosystem offers greater flexibility, features and speed for these kinds of problems. The **DenseGillespieAlgorithm.jl** is optimised for high-dimensional systems. For further insight into the simulation framework and its applications, please refer to the accompanying [documentation](https://roccminton.github.io/DenseGillespieAlgorithm.jl/stable).

### Example: SIR 

```julia
using DenseGillespieAlgorithm

# Define the initial state of the system (e.g., number of molecules)
initial_state = [100, 50]

# Define the reactions (reaction rates and species interactions)
reactions = [
    (rate=0.1, reactants=[1, 0], products=[0, 1]),  # Reaction 1
    (rate=0.05, reactants=[0, 1], products=[1, 0])  # Reaction 2
]

# Simulation parameters
tmax = 100.0  # Maximum simulation time

# Run the simulation
result = simulate_gillespie(initial_state, reactions, tmax)

# Analyze or plot the result (example with a simple print)
println(result)
```
### Key Function
`run_gillespie`: This is the main function that runs the simulation, given the initial state, the affect and the rates function and the maximum time to simulate. The arguments are as follows:
  - `time`: Time intervall for simulation, eg. a range like 1:100
  - `n_0`: The initial population state. This gets overwritten during the simulation
  - `par`: Additional parameters. This is handed to the execute and rates functions and gives extra flexibility
  - `execute!`: Execute function. This gets and index i that determines the event, the population state and the parameter `par` as input.
  - `rates!`: Calculates and overwrites the `initrates`. Therefore gets the rates, the population state and the parameter `par` as input.
  -  `initrates`: The inital rates. 
  -  `population_history`: Empty population history, to save the simulation.
  -  `hstart`: time shifts simulation for parameter changes _(optional, default = 0)_
  -  `statistic!`: addititonal statistic function if not the whole population state should be safed _(opitonal)_

## Contributing

Contributions are welcome! If you encounter any issues or would like to suggest features, please feel free to open an issue or submit a pull request.

To contribute:

1. Fork the repository.
2. Create a new branch (`git checkout -b feature/my-new-feature`).
3. Commit your changes (`git commit -am 'Add some feature'`).
4. Push to the branch (`git push origin feature/my-new-feature`).
5. Create a new Pull Request.

## License

This package is open-sourced under the MIT license. See the [LICENSE](https://github.com/roccminton/DenseGillespieAlgorithm.jl/blob/main/LICENSE) file for details.

---

### Additional Notes

- **Dependencies**: All required dependencies will be automatically installed when the package is added via Julia's package manager. Note for the direct import from github, that the package depends on the following packages:
  - Random
  - Distributions
  - SparseArrays
  - ProgressMeter
- **Support**: For support, feel free to open an issue or reach out via the JuliaLang Discourse.

Happy coding with **DenseGillespieAlgorithm.jl**!
