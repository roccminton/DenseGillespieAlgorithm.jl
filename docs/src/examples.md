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

We present a specific case study of an adaptive dynamics model, originally proposed by Dieckmann and Doebeli[^Dieckmann99]. Here the trait space is ``\mathcal{X} = [-1,1] \subset \mathbb{R}``. The birth rate is given by ``b(x)=exp(-x^2/2\sigma^2_b)`` for some ``\sigma_b > 0``. The death rate is constant ``d(x) = d`` and the competion between individuals depends only on theri distance by ``c(x,y)=exp(-(x-y)^2/2\sigma^2_c)`` for some ``\sigma_c > 0``. Moreover the mutation kernel, that chooses the new trait of an offspring at birth is a Gaussian law with mean ``0`` and variance ``0.1`` conditioned to ``[-1,1]``.

[^Dieckmann99]:U. Dieckmann, M. Doebeli, On the origin of species by sympatric speciation. _Nature_ 400:354-357, 1999

The next step is to begin the implementation of this model, starting with the rates function.


```julia
using Distributions

#birth rate
b(x, σ) = exp(-x^2 / (2σ^2))
#death rate
d(x,d) = d
#competiton kernel
c(x, y, σ, K) = inv(K) * exp(-(x - y)^2 / (2σ^2))
#mutation kernel
mutation(x) = rand(truncated(Normal(x, 0.1), -1, 1))

```

Given that we anticipate a relatively limited number of distinct traits to be present at any given time, but a considerable number of representatives for any given trait that we elect to implement this model with, we have opted to utilise dictionaries. Each trait is a key within the dictionary, with the value being a triple consisting of the size of the subpopulation and its intrinsic birth and death rate. By saving the birth and death rate, the need for repeated recalculation of the same rate in each step is avoided; instead, the rate is simply read from the dictionary.
To illustrate,starting in a monomorphic equilibrium at the boudary ``x_0 = -1``, the initial population state would be as follows.

```julia
x0 = -1.0

n0 = Dict(
        x0 => [
                (b(x0, σ_b) - d(x0,d)) / c(x0,x0,σ_c, K),
                b(x0, σ_b),
                d(x0,d)
                ]
        )
```

In this manner, the rate values for each trait are stored in a cache once they are incorporated into the population. A similar approach is employed for the competition rates between individuals, with a dedicated dictionary being established to accommodate the various competition rates. In order to establish the competition dictionary, it is necessary to define the following function.

```julia
#Generate a cach dictionary for all competition rates between individuals from the population state ps
function generatecompdict(ps,competition)
    IndividualType = keytype(ps)
    Individuals = collect(keys(ps))
    #generate empty dictionary
    C = Dict{
        Tuple{IndividualType,IndividualType},
        Real
        }()
    #populate dictionary
    for x in keys(ps), y in keys(ps)
        C[(x,y)] = competition(x,y)
    end
    return C
end
```
The cached values will be passed to the functions via the parameter variable. Additionally, the birth, death, mutation and competition functions, with their fixed parameter values, will be stored there. Furthermore, it is necessary to adjust all model parameters, including the variances of the Gaussian birth and competition rates, the population size, and the time frame.
```julia
t = 0:1000

par = (
        birth = x -> b(x, 0.9),
        death = x -> d(x,0.0),
        competition = (x, y) -> c(x, y, 0.8, 1000),
        mutate = mutation,
        μ = 0.00015,
        K = 1000,
        compdict = generatecompdict(n0,(x, y) -> c(x, y, 0.8, 1000)),
        historylength = length(t)
)
```
!!! info "History length"
    As the population history will be stored in the dictionary, it is necessary to inform the algorithm of the duration of the simulation. To this end, the field "historylength" must be added to the parameter variable.

The next step is to define the rates function. In this case, the rates are also provided as a dictionary. Each subpopulation has two rates: a birth rate and a death rate. These are calculated from the cache and written to the dictionary.

```julia
#define rates function
function rates!(rates::Dict,ps::Dict,pr)
    #iterate through current population
    for (x,vₓ) in ps
        #size of subpopulation
        nₓ = vₓ[1]
        #check if rates are already cached, if not do so
        !haskey(rates,x) && (rates[x] = valtype(rates)(undef,2))
        #birthrate n_x * b(x)
        rates[x][1] = nₓ*vₓ[2]
        #deathrate n_x * (d(x) + Σ c(x,y) n_y)
        rates[x][2] = nₓ* vₓ[3]
        for (traittuple,c) in pr.compdict
            t₁,t₂ = traittuple
            t₁ == x && (rates[x][2] += nₓ * ps[t₂][1] * c)
        end
    end
end
```
The process of adding a new trait to the population at birth is rendered challenging by the presence of extensive caching, particularly in relation to competition rates. Consequently, a preliminary function is first devised to facilitate the addition of new traits to the population, prior to the implementation of the `birth!` and `death!` functions.

```julia
#add a new trait to current population
function addnewtrait!(ps,rates,par,trait)
    #add to population state
    ps[trait] = [par.diff,par.birth(trait),par.death(trait)]
    #set competition
    for other_trait in keys(ps)
        par.compdict[(trait,other_trait)] = par.competition(trait,other_trait)
        par.compdict[(other_trait,trait)] = par.competition(other_trait,trait)
    end
end
```
The `birth!` and `death!` functions can now be defined with relative ease and combined into a single `execute!` function.
```julia
function birth!(ps, rates, par, trait)
    #Birth with or without mutation
    if par.μ > 0.0 && rand() ≤ par.μ
        #mutate to new type/species and add to species
        new_trait = par.mutate(trait)
        #setup the size of the new type
        if haskey(ps,new_trait)
            ps[new_trait][1] += par.diff
        else
            addnewtrait!(ps,rates,par,new_trait)
        end
    else
        ps[trait][1] += par.diff
    end
    nothing
end

function death!(ps,trait,pr)
    ps[trait][1] -= pr.diff
end

function execute!(i,trait,ps,rates,pr)
    if i==1
        birth!(ps, rates, pr, trait)
    elseif i==2
        death!(ps,trait,pr)
    else
        error("Index Error: Unknown event #$i")
    end
end
```
To initiate the simulation, it is merely necessary to establish an empty rates dictionary and population history, and then to execute the [`run_gillespie!`](@ref) function.
```julia
#empty population history
hist = Dict(x=>zeros(eltype(valtype(n0)),length(t)) for x in keys(n0))
#empty rates dictionary (gets populated in first iteration)
initrates = Dict{keytype(n0),Vector{Real}}()

#execute the simulation
run_gillespie!(
        t,
        n0,
        par,
        execute!,
        rates!,
        initrates,
        hist
)
```
To observe the findings, the size of the subpopulations is plotted over time, with the different traits represented by varying colours.
```julia
using Plots

#setup plot
p=plot(legend=false)

#function to determine the color of the trait
function c(x)
    #find the biggest and smallest key in the population history
    min = min(keys(hist)...)
    max = max(keys(hist)...)
    
    #if there ever has been only one trait return 1 otherwise a color inbetween 
    if min == max
        return 1
    else
        return floor(Integer,((x-min)/(max-min)) * length(cgrad(:thermal))-1) + 1
    end
end

for (x,his_x) in history
    plot!(p,time,his_x,color=cgrad(:thermal).colors.colors[c(x)])
end

p
```
![Simulation result with a mutation rate of 1/K](https://roccminton.github.io/images/TSS.png)

!!! warning "Mutation rate"
    In this scenario, the runtime of the algorithm is highly dependent on the mutation rate. An increase in the mutation rate results in a greater number of different traits. This implementation with dictionaries is most suited for a small number of traits being alive at the same time. However, if the mutation rate is increased to levels of frequent mutation, it is recommended that dictionaries are not used, but instead vectors should be employed for saving the data. The following example demonstrates this technique.

!!!  tip "Population size"
    Nevertheless, increasing the population size in this scenario does not significantly prolong the runtime of the algorithm. This is an advantage of using dictionaries and caching the competition. However, this approach is only effective when the mutation rate is scaled with the population size (as demonstrated in the above example).

!!! info "Empty cache"
    It should be noted that the algorithm performs regular checks for subpopulations in the dictionary with a population size of zero. In the event that such subpopulations are identified, they are removed in order to prevent an excessive expansion of the dictionary. This process is carried out by the [`DenseGillespieAlgorithm.dropzeros!`](@ref) function.

## High-dimensional model
The final example we will present is the most complex. We implement a model to analyse the dynamics of complete recessive lethal diseases. Each disease is triggered by the mutation of a gene and is expressed only in a homozygous state. Therefore, the traitspace for this model is ``\mathcal{X}=\{0,1\}^{2\times N}`` where ``N`` is the number of genes. A detailed description and results of numerous simulations with this exact framework can be found [here](https://arxiv.org/abs/2406.09094)[^LaRocca24].

[^LaRocca24]:L. A. La Rocca, K. Gerischer, A. Bovier, and P. M. Krawitz. Refining the drift barrier hypothesis: a role of recessive gene count and an inhomogeneous muller's ratchet, 2024

Individuals expressing a disease are excluded from the mating process. At birth, each individual randomly selects a fit partner from the population. Following the process of recombination, whereby the diploid genetic information is reduced to a haploid zygote incorporating crossover events, the gametes of the two parents fuse to form a new offspring. New mutations emerge at a constant rate. 

Given that there are ``2^{2N}`` potential configurations with interactions between them, it is not feasible to enumerate them all prior to the start of the simulation. 
Furthermore, it is of no particular interest to ascertain the precise genetic configuration of the entire population. Typically, one is only concerned with summary statistics, such as the mutation burden (the average number of mutations per individual) and the prevalence (the fraction of individuals affected by a disease).
It is therefore only these statistics that will be retained for subsequent analysis. However, for the propagation of the population dynamics, it is essential to have access to the exact configurations.
To be more precise, the total birth and death rates can be calculated via the summary statistics, which we utilise. However, in order to employ an offspring, the configurations are required.


