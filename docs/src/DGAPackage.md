# The DenseGillespieAlgorithm

This package implements a version of the Gillespies algorithm that performs exact stochastic simulations for dense problems.
The Gillespie algorithm[^Gillespie76], introduced by Daniel Gillespie in 1976, is a fundamental tool for simulating the time evolution of systems with discrete, stochastic events, particularly in contexts like biochemical reactions and population dynamics. Its applications are particularly prevalent in contexts such as biochemical reactions and population dynamics. The Gillespie Algorithm is employed to simulate the behaviour of systems wherein reactions or events occur at random intervals. The algorithm generates a sequence of events and their timings by first calculating the rates at which different events or reactions occur. Subsequently, the time until the next event is determined based on these rates, and the type of event that occurs next is selected according to its probability. In the final step, the system state is updated based on the event, and the process is repeated.

[^Gillespie76]: D.T. Gillespie. A general method for numerically simulating the stochastic time evolution of coupled chemical reactions. _Journal of Computational Physics_, 22(4):403-434, 1976

The Gillespie algorithm is a highly renowned and widely utilised technique across diverse communities and ecosystems. A particularly efficient, flexible and comprehensive implementation can be found in the [JumpProcess.jl](https://docs.sciml.ai/JumpProcesses/stable/) package within the SciML ecosystem. We strongly recommend the use of this framework wherever feasible. 

However, the majority of implementations of the Gillespie algorithm require prior knowledge of all potential types and all reactions between these types before the reaction commences. A classic illustration of this is the SIR model (see [Examples](@ref)).
The objective of our implementation in this package is to eliminate this restriction and permit the consideration of both high-dimensional systems, where the precise interactions between every conceivable combination are theoretically possible but practically infeasible, and additionally, systems where the trait space is uncountable, such as the real line. In both cases, the number of distinct traits that are present at any given time is finite, given that the population size is limited. However, new types emerge during the course of the simulation, and the interactions between these types are determined by their specific characteristics. 

## Manual Outline
```@contents
Pages = ["manual.md","examples.md","perform.md","index.md"]
```

## Index

```@index
Modules = [DenseGillespieAlgorithm]
Order   = [:constant, :type, :function, :macro]
```