# Manual

The DenseGillespieAlgorithm framework is designed to assist researchers in simulating their complex models in an exact stochastic manner. It is the responsibility of the user to implement all model-specific functions, such as those pertaining to birth and death events or rate functions. Once this has been done, the framework executes the Gillespie Algorithm and saves the population history. The following section provides an overview of the main function of this package.

## Installation

### Install from GitHub

You can install the package directly from this GitHub repository:

```julia
using Pkg
Pkg.add("https://github.com/roccminton/DenseGillespieAlgorithm.jl")
```

### Install from Julia

Once the package is registered in the official Julia package registry, you can install it via:

```julia
using Pkg
Pkg.add("DenseGillespieAlgorithm")
```

This will install the latest stable version of the package and all required dependencies.

!!! note "Package dependencies"
    When loading the package directly from GitHub, the following packages must be available `Random, Distributions, ProgressMeter, SparseArrays`

## Setting up the model functions

The initial step is to define all interaction functions for the model. In population models, these are typically limited to two: birth and death. However, there is no upper limit on the number of interactions that can be included. As the number of fundamentally different interactions increases, the efficiency of the algorithm is reduced. The framework is speciallised to a small number of different events.

Subsequently, all the interaction functions should be incorporated into a single execute function. This function must accept three inputs: an index, the current population state, and the model parameter. The index specifies which of the defined events should be executed. The current population state is then modified by the event functions. 

```julia
    execute!(i,x,par)
```

Next we need to define the rates function. There must be as many rates as there are interaction events. Therefore the variable `initrates` is usually a `Vector` with as many entries as there are events. The rates function takes as input the initial rates, the current population state and the additional model parameters. The function should calculate the rates according to the population state and modify the `initrates` accordingly.

```julia
    rates!(initrates,x,par)
```

!!! note "Function names"
    The function name may be designated as desired, as they are passed to the core function. The nomenclature is inconsequential.

!!! danger "Function signature"
    Nevertheless, it is crucial to maintain the original function signature, which entails retaining the sequence and the number of arguments as they are called within the algorithmic structure.

!!! tip "Parameter variable"
    There are no restrictions on the parameter variable `par`. Any additional information used to calculate rates and to change the current population state can be added to the parameter element that is passed through all functions. For example, if you want to know the current time of the simulation within the functions you run for time-inhomogeneous models, you could add this to your model parameter.

## Setting up the model parameter, population history and initial population
The final step before running the simulation is to define the model parameters, including the time horizon of the simulation and the initial population state, as well as a blank population history.

The time horizon of the simulation is tipically a `UnitRange`, but can be anything that can be enumerated. 
The type of initial population state should correspond to the functionalities defiend in the function `rates!` and `execute!` as they use and modify this type. 
The empty population history should also match the type of population state, as it will be copied into the population history. In addition, if the population history is a vector or matrix, it should be at least as long as the time horizon.
You can customise the saving process with your own Statistics! function. In this case, you will have to adapt the coupling history to the functionalities of this function. For more details, see [Customized statistics](@ref)

## Execute the simulation
With everything in place, it is time to run the simulations. To do this, call the `run_gillespie!` function from the package. 

```@docs
run_gillespie!
```

Once the simulation has reached its conclusion, the modified population history is returned for further analysis.

## Customized statistics
For many high-dimensional models, the exact configuration at any given time is too much information. In many cases only summery statistics are needed. To avoid accumulating too much data during the runtime of the algorithm that is not needed afterwards, you can define your own `statistics!` function. In this case, only the information you want to collect is stored for further analysis.

As for the `rates!` and `execute!` functions, the function signature is of particular significance. The function accepts as input the population history, which is modified by the function and the current time index, hence the index at which the statistics of the current population state are saved. Additionally, the current state and the model parameter are required.

```julia
    statistics!(population_history,t,x,par)
```

