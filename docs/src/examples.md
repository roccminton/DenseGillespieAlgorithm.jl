# Examples

Three illustrative examples are provided. The first is a minimal working example, designed to facilitate the initial implementation of the framework on the user's machine. The second is a slightly more advanced example, which illustrates the use of an uncountable trait space. The third is a highly complex example, which demonstrates the comprehensive versatility of the package.

## SIR-Model

The SIR model is a tree-dimensional model that is used to model infectious diseases. It is a simple model that assumes that individuals can be placed into one of three categories: susceptible, infected, or recovered. Infected individuals can infect susceptible individuals through random interactions. After becoming infected, individuals can recover and become immune. For more details, see for example [here](https://people.wku.edu/lily.popova.zhuhadar/).

The SIR (susceptible-infected-recovered) model is
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

# Define the reactions (reaction rates and species interactions)
function rates!(rates,x,par)
    #rate of infection
    rates[1] = par.β * x[1] * x[2]
    #rate of recovery
    rates[2] = par.γ * x[2]
    nothing
end

#------
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

## High-dimensional model
