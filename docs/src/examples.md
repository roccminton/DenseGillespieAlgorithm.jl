# Examples

Three illustrative examples are provided. The first is a minimal working example, designed to facilitate the initial implementation of the framework on the user's machine. The second is a slightly more advanced example, which illustrates the use of an uncountable trait space and caching for performance improvement in a particular use case. The third is a highly complex example, which demonstrates the comprehensive versatility of the package. 

## 1. SIR-Model

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

## 2. Continuous trait space
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
function rates!(rates::Dict,ps::Dict,par)
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
        for (traittuple,c) in par.compdict
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
    ps[trait][1] -= par.diff
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

#setup plot
p=plot(legend=false)

for (x,his_x) in history
    plot!(p,time,his_x,color=cgrad(:thermal).colors.colors[c(x)])
end

p
```
![Simulation result with a mutation rate of 1/K](https://roccminton.github.io/images/TSS.png)

Simulation results with a mutation rate of ``1/K`` where ``K`` is the carrying capacity.

!!! warning "Mutation rate"
    In this scenario, the runtime of the algorithm is highly dependent on the mutation rate. An increase in the mutation rate results in a greater number of different traits. This implementation with dictionaries is most suited for a small number of traits being alive at the same time. However, if the mutation rate is increased to levels of frequent mutation, it is recommended that dictionaries are not used, but instead vectors should be employed for saving the data. The following example demonstrates this technique.

!!!  tip "Population size"
    Nevertheless, increasing the population size in this scenario does not significantly prolong the runtime of the algorithm. This is an advantage of using dictionaries and caching the competition. However, this approach is only effective when the mutation rate is scaled with the population size (as demonstrated in the above example).

!!! info "Empty cache"
    It should be noted that the algorithm performs regular checks for subpopulations in the dictionary with a population size of zero. In the event that such subpopulations are identified, they are removed in order to prevent an excessive expansion of the dictionary. This process is carried out by the [`DenseGillespieAlgorithm.dropzeros!`](@ref) function.

!!! info "Switch-off caching"
    It is possible to implement the same dynamics without the caching of competition rates. In this case, the only necessary modification is to alter the `for`-loop in the `rates!` function, which iterates over all pairs of tuples to
    ```julia
    for t₁ in keys(ps)
        rates[x][2] += nₓ * ps[t₁][1] * par.competition(x,t₁)
    end
    ```
    and delete the `for`-loop over all other traits in the `addnewtrait!` function. In instances where the number of distinct traits is considerable, this approach is advised.

## 3. High-dimensional model
The final example we will present is the most complex. We implement a model to analyse the dynamics of complete recessive lethal diseases. Each disease is triggered by the mutation of a gene and is expressed only in a homozygous state. Therefore, the traitspace for this model is ``\mathcal{X}=\{0,1\}^{2\times N}`` where ``N`` is the number of genes. A detailed description and results of numerous simulations with this exact framework can be found [here](https://arxiv.org/abs/2406.09094)[^LaRocca24].

[^LaRocca24]:L. A. La Rocca, K. Gerischer, A. Bovier, and P. M. Krawitz. Refining the drift barrier hypothesis: a role of recessive gene count and an inhomogeneous muller's ratchet, 2024

Individuals expressing a disease are excluded from the mating process. At birth, each individual randomly selects a fit partner from the population. Following the process of recombination, whereby the diploid genetic information is reduced to a haploid zygote incorporating crossover events, the gametes of the two parents fuse to form a new offspring. New mutations emerge at a constant rate. 
It is assumed that the intrinsic death rate and the competition among all individuals are identical. Only the birth rate is reduced to zero for infected individuals.

Given that there are ``2^{2N}`` potential configurations with interactions between them, it is not feasible to enumerate them all prior to the start of the simulation. 
Furthermore, it is of no particular interest to ascertain the precise genetic configuration of the entire population. Typically, one is only concerned with summary statistics, such as the mutation burden (the average number of mutations per individual) and the prevalence (the fraction of individuals affected by a disease).
It is therefore only these statistics that will be retained for subsequent analysis. However, for the propagation of the population dynamics, it is essential to have access to the exact configurations.
To be more precise, the total birth and death rates can be calculated via the summary statistics, which we utilise. However, in order to employ an offspring, the configurations are required.

In order to enhance the efficacy of the algorithm, it is essential to make extensive use of the parameter variable, which is passed to all relevant functions. The intrinsic configuration is stored therein in an optimal way, while the population state encompasses only the summary statistics and the total population size, as these are the necessary components for computing the event rates.

The initial stage of the process entails the establishment of all model-specific parameters, including the constant birth, death, and competition rates; the expected number of mutations at birth, denoted by μ; the number of recessive genes; the initial population size; and the recombination rate.

```julia
using SparseArrays
using Random

t = 0:1000

par = (
    birth = 1.0,
    death = 0.9,
    competition = 0.1 / 1000,
    μ = 0.1,
    Nloci = 100,
    K = 1000,
    recombination = 0.01
)

n0 = Dict(
    "PopSize" => par.K,
    "Ill" => 0,
    "ML" => 0
)
```

The only two possible events are birth and death, the `rates!` function can be expressed as follows: 

```julia
function rates!(rates, ps, par)
    #linear birth for all propagable individuals
    rates[1] = par.birth * (ps["PopSize"]-ps["Ill"])
    #uniform logistic death
    rates[2] = ps["PopSize"] * par.death + ps["PopSize"] * (ps["PopSize"] - 1) * par.competition
    nothing
end
```
The definition of the function that executes the birth is a more challenging undertaking, given that it involves three mechanisms: mating, recombination and mutation. 
Nevertheless, prior to an explication of the implementation of the execute! function, it is necessary to describe the means by which the population configuration is saved internally.
Since the size of any given trait is relatively large (``2N`` bytes), and since we anticipate a significant number of different traits, but a limited population size, we have chosen to construct a vector of traits that is as large as the expected population size, with additional space for fluctuations. This vector will store all the traits. Furthermore, a dictionary of indices is maintained, which points to the indices of traits in the vector. The dictionary differentiates between traits that are either alive and healthy, or alive but ill (thus expressing the disease and unable to reproduce), or that are not part of the current population. The aforementioned free traits can then be modified if new offsprings are born, eliminating the necessity of initiating a new ``2 \times N`` matrix of `Bool`s each time. This method allows for the saving of a considerable amount of memory.
The production of new traits is dependent upon the absence of free indices. Upon the death of a fey individual, the index is released into the group of free indices, where it may be reborn as a new trait at a future point in time.
In order to initiate the trait vector, which encompasses all individual genetic configurations, we have implemented a function that takes the initial population state as an input and draws a possible trait configuration from it. 
```julia
#empty genetic configuration
emptytraits(Nloci,T=Bool) = [spzeros(T,Nloci),spzeros(T,Nloci)]

#produce trait collection from population state
function inittraits(par,n0)
    #Setup healthy genetic information
    locs = 1:par.Nloci
    #Generate healty population with some buffer for fluctuations
    traits = [emptytraits(par.Nloci) for _ in 1:round(Int,par.K + sqrt(par.K))]
    #add two mutations to completely healthy individuals to get the required number of ill individuals
    for i in 1:n0["Ill"]
        l = rand(locs)
        traits[i][1][l] = 1
        traits[i][2][l] = 1
    end
    individuals = 1:n0["PopSize"]
    #add the remaining mutaions to the population to get the required mutation load
    for i in n0["Ill"]+1:n0["ML"]-2*n0["Ill"]
        #choose random individual and location
        ind = rand(individuals)
        l = rand(locs)
        #recoose random individual and location if the individual has already a mutation
        #at that locationo or at the homologe gene
        while traits[ind][1][l]+traits[ind][2][l] ≠ 0
            ind = rand(individuals)
            l = rand(locs)
        end
        traits[ind][rand(par.choosecopyfrom)][l] = 1
    end
    return traits
end
```
Subsequently, both the traits and the corresponding index dictionary are incorporated into the parameter variable. Additionally, other necessary elements are included, allowing for their reuse rather than generation on each occasion. These include a vector of random numbers for mate selection, the mutation distribution, the distribution of mutation locations and a unit range for gene segment selection.
```julia
par = (par...,
    rndm = Vector{Int}(undef,2),
    MutationsPerBirth = Poisson(par.μ),
    MutationLocation = 1:par.Nloci,
    traits = inittraits(par,n0),
    indices = Dict(
        "healthy" => collect(n0["Ill"]+1:n0["PopSize"]),
        "ill" => collect(1:n0["Ill"]),
        "free" => collect(n0["PopSize"]+1:round(Int,par.K + sqrt(par.K)))
    ),
    historylength = length(t),
    choosecopyfrom = 1:2,
)
```
The initial step on implementing the `birth!` function is to implement a function that establishes the crossover breakpoints for recombination at random, in accordance with the specified recombination rate. The function initially draws the number of crossover breakpoints from a _Poisson_ distribution, and subsequently selects the position at random from among all ``N-1`` positions. The resulting vector of `UnitRange` segments is then returned.

```julia
#output for recombination rate 1
fullreccuts(par) = [i:i for i in 1:par.Nloci]
#output for recombination rate 0
noreccuts(par) = [1:par.Nloci]

function reccuts(par)
    if par.recombination == 1
        return fullreccuts(par)
    elseif par.recombination == 0
        return noreccuts(par)
    else
        #draw number of chromosome cuts
        ncuts = rand(Poisson(par.recombination*par.Nloci))
        #equals full recombination
        ncuts ≥ par.Nloci - 1 && return fullreccuts(par)
        #equals no recombination
        iszero(ncuts) && return noreccuts(par)
        #otherwise produce individual segments at random
        cutsat = sort!(sample(1:par.Nloci-1,ncuts,replace=false))
        ccuts = [1:cutsat[1]]
        for i in 2:length(cutsat)
            push!(ccuts,cutsat[i-1]+1:cutsat[i])
        end
        push!(ccuts,cutsat[end]+1:par.Nloci)
        return ccuts
    end
end
```
Subsequently, once two parents have been identified for mating, the process of the offspring generation is defined. This encompasses recombination, mating and mutation. The new genetic configuration of the offspring is stored at a designated index.
```julia
function offspring!(offspring_index, par, n_mut)
    #randomly recombine the parental genetic information
    #first for one then for the other parent
    for i in par.choosecopyfrom # =1:2
        #randomly choose one copy for each chromosome/gene block
        ccuts = reccuts(par)
        choosecopy = rand(par.choosecopyfrom,length(ccuts))
        for (r,chromosome) in enumerate(ccuts)
            view(par.traits[offspring_index][i],chromosome) .=
                view(par.traits[par.rndm[i]][choosecopy[r]],chromosome)
        end
    end
    #add n_mut mutations to random positions mutation
    #if there are no mutations to add skip the mutation process
    if n_mut > 0
        for _ in 1:n_mut
            par.traits[offspring_index][rand(par.choosecopyfrom)][rand(par.MutationLocation)] = 1
        end
    end
    nothing
end
```
The final step before integrating all components into a unified system is to implement a function that updates the population state following the generation of the offspring's genetic configuration.It is therefore necessary to ascertain whether the offspring in question exhibits a mutation in a homogeneous state, which would render it unsuitable for reproduction. Furthermore, the mutation burden of the offspring must be calculated.
```julia
#check if configuration has homogeneous mutation
ispropagable(a::Vector,Nloci) = ispropagable(a)
function ispropagable(a::Vector)
    for (i,p) in enumerate(a[1])
        isone(p) && isone(a[2][i]) && return false
    end
    return true
end
#calculate mutation load
mutationload(a::Vector) = sum(sum(svec) for svec in a)
#modify population state at birth of offspring
function updateps_birth!(ps,par,offspring_index)
    if ispropagable(par.traits[offspring_index],par.Nloci)
        push!(par.indices["healthy"],offspring_index)
    else
        ps["Ill"] += 1
        push!(par.indices["ill"],offspring_index)
    end
    ps["PopSize"] += 1
    ps["ML"] += mutationload(par.traits[offspring_index])
end
```
The aforementioned processes are unified in the `birth!` function, which initially identifies potential parents for mating, subsequently generates offspring, and finally updates the population state. 
```julia
#execute the addition of an individual
function birth!(ps, par)
    #choose two genetic configurations to mate
    rand!(par.rndm,par.indices["healthy"])
    #clean up parental configurations
    for i in par.choosecopyfrom, j in par.choosecopyfrom
        dropzeros!(par.traits[par.rndm[i]][j])
    end
    #select free index for offspring
    if isempty(par.indices["free"])
        offspring_index = length(par.traits) + 1
        push!(par.traits,emptytraits(par.Nloci))
    else
        offspring_index = pop!(par.indices["free"])
    end
    #generate offsprings genetic configuration
    offspring!(offspring_index, par, rand(par.MutationsPerBirth))
    #add the individual to the current population state dictionary
    updateps_birth!(ps,par,offspring_index)
end
```
The removal of an individual at death is a relatively straightforward process. It merely entails freeing the index and updating the population state in accordance with the relevant changes.
```julia
#modify population state at death of an individual
function updateps_death!(ps,par,fey_index)
    ps["PopSize"] -= 1
    ps["ML"] -= mutationload(par.traits[fey_index])
end
#execute the removal of an individual
function death!(ps, par)
    #choose fey
    if rand()<=ps["Ill"]/ps["PopSize"]
        fey_index = popat!(par.indices["ill"],rand(1:ps["Ill"]))
        ps["Ill"] -= 1
    else
        fey_index = popat!(par.indices["healthy"],rand(1:(ps["PopSize"]-ps["Ill"])))
    end
    #add fey to gravejard
    push!(par.indices["free"],fey_index)
    #update population state
    updateps_death!(ps,par,fey_index)
end
```
As in the preceding cases, it is necessary to collate both the `birth!` and `death!` functions into a single `execute!` function.
```julia
function execute!(i,ps,par)
    if i == 1
        birth!(ps,par)
    elseif i == 2
        death!(ps,par)
    else
        error("Unknown event index: $i")
    end
end
```
Finally, the simulation can be executed following the initialisation of both a blank rates vector and a blank population history.
```julia
#setup empty rates vector
initrates = Vector{typeof(par.birth)}(undef,2)

#setup empty population history
hist = Dict(x=>zeros(valtype(n0),length(t)) for x in keys(n0))

run_gillespie!(
        t,n0,
        par,
        execute!,
        rates!,
        initrates,
        hist
        )
```
In order to facilitate the analysis of the data, we present a straightforward graphical representation.
```julia
using Plots
#plot the prevalence
plot(t,hist["Ill"] ./ hist["PopSize"],color=:orange,label="")
#plot the population size (without axis ticks)
plot!(twinx(),hist["PopSize"],color=:gray,yticks=false,label="")
#plot the mutation burden
plot!(twinx(),hist["ML"] ./ hist["PopSize"],color=:red,label="")
```
!!! warning "Population extinction"
    In the event of the population becoming extinct, the aforementioned method will result in an error, as the division of zero will occur. To circumvent this issue, it is recommended to first invoke the following function on the data set.
    ```julia
    replace_NaN(v) = map(x -> isnan(x) ? zero(x) : x, v)
    ```
![Simulation result showing the mutation load, prevalence and population size](https://roccminton.github.io/images/MLP.png)

The grey line represents the population size, which is not represented by any axis. On the left y-axis, the prevalence is represented by the yellow line, while the mutation burden is shown by the red line on the right y-axis.

### Custom `statistic!` function
Thus far, the only function employed for the purpose of saving the population history was the built-in [`DenseGillespieAlgorithm.saveonestep!`](@ref) function, which in this case saves the mutation burden, prevalence, and population size over time. However, given the intricate population structure of the model, it may be beneficial to consider saving additional statistics that extend beyond mere numbers over time. For this example, we are interested in saving the allele frequencies of the mutated allele per position over time. This necessitates the definition of a custom `statistic!` function.

The initial step is to incorporate an additional function call into the existing functions `updateps_death!` and `updateps_birth!` of the form `updatestats_death!(ps, par, fey_index)` and `updatestats_birth!(ps, par, offspring_index)`, respectively. This allows us to modify the new statistic that we wish to utilise at each event, rather than recalculating it from the current population state after every full time step, which is the usual process.

!!! tip "Empty functions"
    In the event that there is no intention to update the statistical data at each stage, it is nonetheless recommended that the function call `updatestats_event!` be retained in order to ensure the flexibility and reusability of the code. In the absence of a required function, the implementation of a generic function of the form 
    ```julia
    function updatestats_death! end
    function updatestats_birth! end
    ```
    is sufficient.

In order to enhance the flexibility of the system, a custom data type has been defined to accommodate the various statistical elements associated with the population history. This approach facilitates the incorporation of new statistics, should the need arise.

```julia
#type to hold population history
struct PopHist
    #mutation burden, prevalence and population size
    mlp :: Dict 
    #allele frequencies per postion
    loadpos :: Array
end
```
!!! note "Unmutable struct"
    It is important to note that the variable type of the population history has been defined as unmutable, which may appear counterintuitive at first glance. However, upon closer examination, it becomes evident that the elements within the struct are only generated once and then populated with data. Meanwhile, the container itself (array, dict, etc.) remains unchanged. This allows for the use of a faster and lighter unmutable object. Conversely, if there is a need to modify the fields within the struct, it would be necessary to define it as a `mutable struct`.

We choose to save the allele frequencies of the mutated allele at each position as a ``T\times N``  matrix, where ``T`` is the total length of the simulation. Each column of the matrix represents the allele frequencies for a single time step. In fact, we will save the precise number of mutations per gene, leaving the division by the population size to be performed subsequently, once the simulation has been completed.
In order to utilise the enhanced performance afforded by the addition of the updatestats_event function, it is necessary to create a temporary storage location for the current allele frequencies prior to their final saving to the storage medium for subsequent analysis. Once again, the parameter variable that is passed to every significant function is employed for this purpose.
```julia
#add blank current allele frequencies to parameter variable
par = (
    par...,
    cafs = zeros(Int,par.Nloci)
)
```
Given the type configuration that is added or removed from the population, the adjustment of the current allele frequencies is a relatively straightforward process.
```julia
#add or remove one individual from allele frequency vector
function update_allelefreqs!(af,ind,i)
    af .+= i .* ind[1]
    af .+= i .* ind[2]
end
```
Furthermore, the corresponding functions for birth and death can be developed upon this function.
```julia
updatestats_death!(ps,par,index) = update_allelefreqs!(par.cafs,par.traits[index],-1)
updatestats_birth!(ps,par,index) = update_allelefreqs!(par.cafs,par.traits[index],+1)
```
The additional statistics have now been incorporated into the system, and the next stage is to save the data at each time step within the specified time horizon into the `PopHist` type. This process is carried out by the following function. Additionally, the function responsible for saving the supplementary statistical data is merged with the one that stores the mutation burden, prevalence, and population size, which were previously saved in the basic example. This allows the creation of the custom `statistic!` function.

```julia
function saveafs!(allelefreqs,index,ps,par)
        view(allelefreqs,index,:) .= par.cafs
end

function statistic!(pophist::PopHist,index,ps,par)
    #save standard statistic
    DenseGillespieAlgorithm.saveonestep!(pophist.mlp,index,ps,par)
    #save additional statistic
    saveafs!(pophist.allelefreqs,index,ps,par)
   end
```
As previously described, the final stage of the process is to set up the initial rates, the initial population history, and then to execute the run_gillespie! function together with the newly defined `statistic!` function as a keyword argument.

```julia
#setup empty rates vector
initrates = Vector{typeof(par.birth)}(undef,2)

#setup empty population history
hist = PopHist(
    Dict(x=>zeros(valtype(n0),length(t)) for x in keys(n0)),
    zeros(Integer,(length(t),par.Nloci,))
    )

run_gillespie!(
        t,n0,
        par,
        execute!,
        rates!,
        initrates,
        hist,
        statistic! =statistic!
        )
```
To analyse the data, one possible approach would be to construct a small GIF that generates plots of the allele frequencies at each position over time.
```julia
using Plots

#calculate frequencies from absolute numbers of mutations
afs = hist.allelefreqs ./ hist.mlp["PopSize"]
#maximal frequencie for axis limit
ymax = maximum(afs)
#create animation
anim = @animate for i in 0:100
    bar(view(afs,i+1,:),ylim=(0,ymax),label="")
end every 10
#show animation
gif(anim)
```

![Animation of allele frequencies over time](https://roccminton.github.io/images/AFs.gif)

!!! tip "Time intervals"
    In certain instances, the requisite statistic may require a considerable amount of memory space or a significant amount of time to calculate. In such cases, it may be more efficient to save and calculate the statistic not in every time step, but rather only after larger intervals. This can be achieved in two ways.
    First, the time horizon provided to the algorithm can be adjusted to a coarser resolution, for instance, `0:10:1000` instead of `0:1000`, resulting in a step size of 10 rather than 1. It should be noted that in such instances, the events continue to occur at the (potentially very small) event rates, but the saving mechanism is executed at each full time step. In this scenario, however, all the statistics that are generated are saved exclusively at the larger time steps.
    Second, in the event that a specific statistic is particularly resource-intensive, it is possible to implement an `if` condition within the `statistics!` function that will then save the statistic only if the time index meets the specified condition.

!!! tip "Snapshots"
    It is a source of considerable frustration to have invested a significant amount of time in running a comprehensive simulation only to realise, upon completion, that an alternative statistic might also warrant examination. Consequently, it was beneficial on occasion to also take "snapshots" at every couple of generations.
    Therefore, a random sample of the population was selected and all the information for that subpopulation was stored. As the population size was reduced by taking a sample from the population and the time horizon was reduced by taking these snapshots on a coarse time grid, the amount of memory required remained within acceptable limits.