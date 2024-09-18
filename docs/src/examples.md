# Examples

Three illustrative examples are provided. The first is a minimal working example, designed to facilitate the initial implementation of the framework on the user's machine. The second is a slightly more advanced example, which illustrates the use of an uncountable trait space. The third is a highly complex example, which demonstrates the comprehensive versatility of the package.

## SIR-Model

The SIR model is a tree-dimensional model that is used to model infectious diseases. It is a simple model that assumes that individuals can be placed into one of three categories: susceptible, infected, or recovered. Infected individuals can infect susceptible individuals through random interactions. After becoming infected, individuals can recover and become immune. For more details, see for example [here](https://people.wku.edu/lily.popova.zhuhadar/).

The initial step is to implement the fundamental interaction functions. In this scenario, two events are occurring: infection and recovery. The objective is to implement these functions in a manner that modifies the population state, which is represented as a vector with three entries, one for each possible state of an individual.

```julia
using Plots
using DenseGillespieAlgorithm

# Define the reactions
function infection!(x)
    x[1] += -1
    x[2] += 1
    nothing
end

function recovery!(x)
    x[2] += -1
    x[3] += 1
    nothing
end
```

The subsequent step is to combine the aforementioned two functions into a single execute function, which will subsequently be provided to the algorithm.

```julia
#Combine all reactions into one execute! function
function execute!(i,x,par)
    if i == 1
        infection!(x)
    elseif i == 2
        recovery!(x)
    else
        error("Unknown event number i = $i")
    end
    nothing
end
```

Furthermore, define the rate at which the events occur, which also depends on the population state. It should be noted that a modification of an existing variable that holds the current rates is necessary. In this case, as there are two events, namely infection and recovery, the rates variable will be a vector with two entries.

```julia
# Define the reactions (reaction rates and species interactions)
function rates!(rates,x,par)
    #rate of infection
    rates[1] = par.β * x[1] * x[2]
    #rate of recovery
    rates[2] = par.γ * x[2]
    nothing
end
```

Prior to commencing the study, it is essential to define all relevant model parameters. These include the interaction rates, the initial population state and the time horizon for the simulation.
Prior to commencing the simulation, it is essential to define all relevant model parameters. These include the interaction rates, the initial population state and the time horizon for the simulation. Additionally, it is necessary to provide the algorithm with an empty population history, which will be populated with data during runtime.

```julia
#Define the model parameter 
par = (
    β = 0.000005,
    γ = 0.005
    )

# Define the initial state of the system 
x0 = [9999,1,0]

# Define the time horizon for the simulation
t = 0:2000

# Initialize population history
hist = zeros(Int,(length(t),3))
```

At this point in the process, all the necessary components have been put in place, and the task can be handed over to the core function of the package. Once the simulation has been executed, the results are plotted.

```julia
# Run the simulation
run_gillespie!(
        t,x0,par,
        execute!,rates!,
        Vector{Float64}(undef,2),hist
        )

# Analyze or plot the result (example with a simple print)
plot(hist,label=["S" "I" "R"])
```
![SIR Plot](https://roccminton.github.io/images/SIR.png)

!!! note "Note"
    While it is feasible to construct such straightforward examples using the DenseGillespieAlgorithm, this is not the typical application. For relatively simple models, the [JumpProcess.jl](https://docs.sciml.ai/JumpProcesses/stable/) package offers greater flexibility and facilitates the implementation process.

## Continuous trait space
In this example, we consider an individual-based model of adaptive dynamics, wherein the trait space is a subset of the real line.
It is therefore impossible to list all the types and interaction rates between them, as there are uncountably many. It is thus necessary to implement the rates and interactions in a dynamic manner. 

In the context of adaptive dynamics models, individuals are characterised by a specific trait, which in this case is a real number. The mortality and fertility rates of individuals are contingent upon this trait. Moreover, competition among individuals is contingent upon the trait in question. Furthermore, at birth, with a probability of μ, the offspring undergoes a mutation and displays a distinct trait in comparison to its parents. 

For further insight into the subject of adaptive dynamics models, we would direct the reader to the [lecture notes](https://www.dropbox.com/scl/fi/8f0pdhq471unbla5dw5xh/LN_SMLS.pdf?rlkey=b0eodhvonuehueuancxeoqloz&st=eqd3ffwj&dl=0) by Anton Bovier.

We present a specific case study of an adaptive dynamics model, originally proposed by Dieckmann and Doebeli[^Dieckmann99]. 

[^Dieckmann99]:U. Dieckmann, M. Doebeli, On the origin of species by sympatric speciation. _Nature_ 400:354-357, 1999


## High-dimensional model
