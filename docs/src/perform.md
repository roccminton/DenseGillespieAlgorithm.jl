# Performance Tips

One of the key benefits of the Gillespie algorithm is its ability to trace a single, precise stochastic trajectory. Nevertheless, in order to achieve this for each individual event, the rates must be calculated and re-calculated whenever there is a change in the population configuration. This makes the algorithm computationally demanding. There are numerous modifications that can be made in order to enhance performance, such as tau-leaping[^Gillespie01]. However, in this section, our objective is to maintain the precision of the stochastic simulation and to identify potential bottlenecks and strategies for optimising the performance of the simulation in its current form.

[^Gillespie01]: D.T. Gillespie. Approximate accelerated stochastic simulation of chemically reacting systems. _Journal of Chemical Physics_, 115(4):1716-1733, 2001