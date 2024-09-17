var documenterSearchIndex = {"docs":
[{"location":"examples.html#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples.html","page":"Examples","title":"Examples","text":"Three illustrative examples are provided. The first is a minimal working example, designed to facilitate the initial implementation of the framework on the user's machine. The second is a slightly more advanced example, which illustrates the use of an uncountable trait space. The third is a highly complex example, which demonstrates the comprehensive versatility of the package.","category":"page"},{"location":"examples.html#SIR-Model","page":"Examples","title":"SIR-Model","text":"","category":"section"},{"location":"examples.html","page":"Examples","title":"Examples","text":"The SIR model is a tree-dimensional model that is used to model infectious diseases. It is a simple model that assumes that individuals can be placed into one of three categories: susceptible, infected, or recovered. Infected individuals can infect susceptible individuals through random interactions. After becoming infected, individuals can recover and become immune. For more details, see for example here.","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"The SIR (susceptible-infected-recovered) model is","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"using Plots\nusing DenseGillespieAlgorithm\n\n# Define the reactions\nfunction infection!(x)\n    x[1] += -1\n    x[2] += 1\n    nothing\nend\n\nfunction recovery!(x)\n    x[2] += -1\n    x[3] += 1\n    nothing\nend\n\n#Combine all reactions into one execute! function\nfunction execute!(i,x,par)\n    if i == 1\n        infection!(x)\n    elseif i == 2\n        recovery!(x)\n    else\n        error(\"Unknown event number i = $i\")\n    end\n    nothing\nend\n\n# Define the reactions (reaction rates and species interactions)\nfunction rates!(rates,x,par)\n    #rate of infection\n    rates[1] = par.β * x[1] * x[2]\n    #rate of recovery\n    rates[2] = par.γ * x[2]\n    nothing\nend\n\n#------\n#Define the model parameter \npar = (\n    β = 0.000005,\n    γ = 0.005\n    )\n\n# Define the initial state of the system \nx0 = [9999,1,0]\n\n# Define the time horizon for the simulation\nt = 0:2000\n\n# Initialize population history\nhist = zeros(Int,(length(t),3))\n\n# Run the simulation\nrun_gillespie!(\n        t,x0,par,\n        execute!,rates!,\n        Vector{Float64}(undef,2),hist\n        )\n\n# Analyze or plot the result (example with a simple print)\nplot(hist,label=[\"S\" \"I\" \"R\"])","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"(Image: SIR Plot)","category":"page"},{"location":"examples.html","page":"Examples","title":"Examples","text":"note: Note\nWhile it is feasible to construct such straightforward examples using the DenseGillespieAlgorithm, this is not the typical application. For relatively simple models, the JumpProcess.jl package offers greater flexibility and facilitates the implementation process.","category":"page"},{"location":"examples.html#Continuous-trait-space","page":"Examples","title":"Continuous trait space","text":"","category":"section"},{"location":"examples.html#High-dimensional-model","page":"Examples","title":"High-dimensional model","text":"","category":"section"},{"location":"perform.html#Performance-Tips","page":"Performance Tips","title":"Performance Tips","text":"","category":"section"},{"location":"perform.html","page":"Performance Tips","title":"Performance Tips","text":"One of the key benefits of the Gillespie algorithm is its ability to trace a single, precise stochastic trajectory. Nevertheless, in order to achieve this for each individual event, the rates must be calculated and re-calculated whenever there is a change in the population configuration. This makes the algorithm computationally demanding. There are numerous modifications that can be made in order to enhance performance, such as tau-leaping[Gillespie01]. However, in this section, our objective is to maintain the precision of the stochastic simulation and to identify potential bottlenecks and strategies for optimising the performance of the simulation in its current form.","category":"page"},{"location":"perform.html","page":"Performance Tips","title":"Performance Tips","text":"[Gillespie01]: D.T. Gillespie. Approximate accelerated stochastic simulation of chemically reacting systems. Journal of Chemical Physics, 115(4):1716-1733, 2001","category":"page"},{"location":"manual.html#Manual","page":"Manual","title":"Manual","text":"","category":"section"},{"location":"manual.html","page":"Manual","title":"Manual","text":"The DenseGillespieAlgorithm framework is designed to assist researchers in simulating their complex models in an exact stochastic manner. It is the responsibility of the user to implement all model-specific functions, such as those pertaining to birth and death events or rate functions. Once this has been done, the framework executes the Gillespie Algorithm and saves the population history. The following section provides an overview of the main function of this package.","category":"page"},{"location":"manual.html#Installation","page":"Manual","title":"Installation","text":"","category":"section"},{"location":"manual.html#Install-from-GitHub","page":"Manual","title":"Install from GitHub","text":"","category":"section"},{"location":"manual.html","page":"Manual","title":"Manual","text":"You can install the package directly from this GitHub repository:","category":"page"},{"location":"manual.html","page":"Manual","title":"Manual","text":"using Pkg\nPkg.add(\"https://github.com/roccminton/DenseGillespieAlgorithm.jl\")","category":"page"},{"location":"manual.html#Install-from-Julia","page":"Manual","title":"Install from Julia","text":"","category":"section"},{"location":"manual.html","page":"Manual","title":"Manual","text":"Once the package is registered in the official Julia package registry, you can install it via:","category":"page"},{"location":"manual.html","page":"Manual","title":"Manual","text":"using Pkg\nPkg.add(\"DenseGillespieAlgorithm\")","category":"page"},{"location":"manual.html","page":"Manual","title":"Manual","text":"This will install the latest stable version of the package and all required dependencies.","category":"page"},{"location":"manual.html","page":"Manual","title":"Manual","text":"note: Package dependencies\nWhen loading the package directly from GitHub, the following packages must be available Random, Distributions, ProgressMeter, SparseArrays","category":"page"},{"location":"manual.html#Setting-up-the-model-functions","page":"Manual","title":"Setting up the model functions","text":"","category":"section"},{"location":"manual.html","page":"Manual","title":"Manual","text":"The initial step is to define all interaction functions for the model. In population models, these are typically limited to two: birth and death. However, there is no upper limit on the number of interactions that can be included. As the number of fundamentally different interactions increases, the efficiency of the algorithm is reduced. The framework is speciallised to a small number of different events.","category":"page"},{"location":"manual.html","page":"Manual","title":"Manual","text":"Subsequently, all the interaction functions should be incorporated into a single execute function. This function must accept three inputs: an index, the current population state, and the model parameter. The index specifies which of the defined events should be executed. The current population state is then modified by the event functions. ","category":"page"},{"location":"manual.html","page":"Manual","title":"Manual","text":"    execute!(i,x,par)","category":"page"},{"location":"manual.html","page":"Manual","title":"Manual","text":"Next we need to define the rates function. There must be as many rates as there are interaction events. Therefore the variable initrates is usually a Vector with as many entries as there are events. The rates function takes as input the inital rates, the current population state and the additional model parameters. The function should calculate the rates according to the population state and modify the initrates accordingly.","category":"page"},{"location":"manual.html","page":"Manual","title":"Manual","text":"    rates!(initrates,x,par)","category":"page"},{"location":"manual.html","page":"Manual","title":"Manual","text":"note: Function names\nThe function name may be designated as desired, as they are passed to the ` function. The nomenclature is inconsequential.","category":"page"},{"location":"manual.html","page":"Manual","title":"Manual","text":"danger: Function signature\nNevertheless, it is crucial to maintain the original function signature, which entails retaining the sequence and the number of arguments as they are called within the algorithmic structure.","category":"page"},{"location":"manual.html","page":"Manual","title":"Manual","text":"tip: Parameter variable\nThere are no restrictions on the parameter variable par. Any additional information used to calculate rates and to change the current population state can be added to the parameter element that is passed through all functions. For example, if you want to know the current time of the simulation within the functions you run for time-inhomogeneous models, you could add this to your model parameter.","category":"page"},{"location":"manual.html#Setting-up-the-model-parameter,-population-history-and-initial-population","page":"Manual","title":"Setting up the model parameter, population history and initial population","text":"","category":"section"},{"location":"manual.html","page":"Manual","title":"Manual","text":"The final step before running the simulation is to define the model parameters, including the time horizon of the simulation and the initial population state, as well as a blank population history.","category":"page"},{"location":"manual.html","page":"Manual","title":"Manual","text":"The time horizon of the simulation is tipically a UnitRange, but can be anything that can be enumerated.  The type of initial population state should correspond to the functionalities defiend in the function rates! and execute! as they use and modify this type.  The empty population history should also match the type of population state, as it will be copied into the population history. In addition, if the population history is a vector or matrix, it should be at least as long as the time horizon. You can customise the saving process with your own Statistics! function. In this case, you will have to adapt the coupling history to the functionalities of this function. For more details, see Customized Statistics","category":"page"},{"location":"manual.html#Execute-the-simulation","page":"Manual","title":"Execute the simulation","text":"","category":"section"},{"location":"manual.html","page":"Manual","title":"Manual","text":"With everything in place, it is time to run the simulations. To do this, call the run_gillespie! function from the package. ","category":"page"},{"location":"manual.html","page":"Manual","title":"Manual","text":"run_gillespie!","category":"page"},{"location":"manual.html#DenseGillespieAlgorithm.run_gillespie!","page":"Manual","title":"DenseGillespieAlgorithm.run_gillespie!","text":"    run_gillespie!(time,n₀,par,execute!,rates!,initrates,population_history[,hstart=0,statistic!])\n\nRun a exact stochastic simulation, return and fill the population_history.\n\nArguments\n\ntime::AbstracVector: time interval for the simulation\nn₀: initial population state\npar: additional parameter (gets passed to execute! and rates!)\nexecute!: execute function\nrates!: rates function\ninitrates: initial rates\npopulation_history: empty population history\nhstart=0: time shift for parameter change (opitonal)\nstatistic!: additional statistic function (optional)\n\nExtended help\n\nNote that n₀,initrates,population_history all three get modified during the simulation\nThe algorithm expects the execute! function to have the following signature   julia   execute!(i::Number,n₀,par)   where the i is the event that gets executed and the population state n₀ gets modified accordingly.   The only exception is when the initrates are given as a dictionary. In that case the signature is execute!(i,trait,n₀,initrates,par), where 'trait'   is the key that is modified.\nThe algorithm expects the rates! function to have the following signature   julia   rates!(initrates,n₀,par)   where the rates get modified according to the current population state given in n₀.\nThe algorithm expects the statistic! function to have the following signature   julia   statistic!(population_hist,t,n₀,par)   where the population history gets modified at position t with the current population state n₀.\nNote that the population_history needs to be accessable via index from 1 to length(time), or if hstart is given from 1+hstart to length(time)+hstart. Unless a specified statistic! function is given.\nNote that the initial population state n₀ must match the population_history in the sense that population_history :: Vector{typeof(n₀)}.  Unless a specified statistic! function is given.\nThe parameter variable par is passed through all functions (execute!,rates!,statistics!), thereby affording the user additional flexibility.\n\n\n\n\n\n","category":"function"},{"location":"manual.html","page":"Manual","title":"Manual","text":"Once the simulation has reached its conclusion, the modified population history is returned for further analysis.","category":"page"},{"location":"manual.html#Customized-Statistics","page":"Manual","title":"Customized Statistics","text":"","category":"section"},{"location":"manual.html","page":"Manual","title":"Manual","text":"For many high-dimensional models, the exact configuration at any given time is too much information. In many cases only summery statistics are needed. To avoid accumulating too much data during the runtime of the algorithm that is not needed afterwards, you can define your own statistics! function. In this case, only the information you want to collect is stored for further analysis.","category":"page"},{"location":"manual.html","page":"Manual","title":"Manual","text":"As for the rates! and execute! functions, the function signature is of particular significance. The function accepts as input the population history, which is modified by the function and the current time index, hence the index at which the statistics of the current population state are saved. Additionally, the current state and the model parameter are required.","category":"page"},{"location":"manual.html","page":"Manual","title":"Manual","text":"    statistics!(population_history,t,x,par)","category":"page"},{"location":"index.html#Public-API","page":"Public API","title":"Public API","text":"","category":"section"},{"location":"index.html","page":"Public API","title":"Public API","text":"Puplic documentation of all internal functions. ","category":"page"},{"location":"index.html#Detailed-API","page":"Public API","title":"Detailed API","text":"","category":"section"},{"location":"index.html","page":"Public API","title":"Public API","text":"Modules = [DenseGillespieAlgorithm]\nOrder   = [:constant, :type, :function, :macro]","category":"page"},{"location":"index.html#DenseGillespieAlgorithm.chooseevent-Tuple{Any, Any}","page":"Public API","title":"DenseGillespieAlgorithm.chooseevent","text":"chooseevent(rates::Vector{Float64},total_rate::Float64)\n\nChoose from the vector of total rates at random one of the indices of the vector according to their rates.\nThe value 0 is returned if the total rates are positive, but too smale to let the evolution continue.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.chooseevent-Tuple{Dict, Any}","page":"Public API","title":"DenseGillespieAlgorithm.chooseevent","text":"chooseevent(rates::Dict,total_rate::Float64)\n\nChoose from the dictionary of total rates at random one of the keys of the dictionary according to their values.\nThe value 0 is returned if the total rates are positive, but too smale to let the evolution continue.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.dropzeros!-Tuple{Any}","page":"Public API","title":"DenseGillespieAlgorithm.dropzeros!","text":"dropzeros!(ps)\n\nDo nothing for non-dictionary inputs.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.dropzeros!-Tuple{Dict{<:Any, <:Vector}}","page":"Public API","title":"DenseGillespieAlgorithm.dropzeros!","text":"dropzeros!(ps::Dict{Any,Vector})\n\nEliminates all key value pairs for which the the firts entry of the vector of the value is zero.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.historylength-Tuple{Vector, Any}","page":"Public API","title":"DenseGillespieAlgorithm.historylength","text":"historylength(population_hisotry,par)\n\nReturn the simulation time based on the length of the population history.\nIf the population history is neither a `Vector` nor a `Matrix` it is assumed that the Parameter has a field called `historylength` that is then returned.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.mainiteration!-Union{Tuple{F3}, Tuple{F2}, Tuple{F1}, Tuple{Any, Any, Any, Any, Any, Any, F1, F2, F3, Any}} where {F1, F2, F3}","page":"Public API","title":"DenseGillespieAlgorithm.mainiteration!","text":"mainiteration!(pop_hist,rates,n0,ct,time,par,ex!::F1,r!::F2,stat!::F3,hstart)\n\nMainiteration of the GillespieAlgorithm for complex models.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.mainiteration!-Union{Tuple{F3}, Tuple{F2}, Tuple{F1}, Tuple{Any, Any, Real, Any, Any, Any, F1, F2, F3, Any}} where {F1, F2, F3}","page":"Public API","title":"DenseGillespieAlgorithm.mainiteration!","text":"mainiteration!(pop_hist,rates,n0::Real,ct,time,par,ex!::F1,r!::F2,stat!::F3,hstart)\n\nMainiteration of the GillespieAlgorithm for OneType Models where the population state is a number and the population history is a vector.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.nexteventandtime-Tuple{Any}","page":"Public API","title":"DenseGillespieAlgorithm.nexteventandtime","text":"nexteventandtime(rates::Vector{Float64})\n\nSample a exponential distributed random variable to determine the time for the next event and calls `choose_event`.\nThe return value is a tuple consiting of the envent index returned by `choose_event` and the time to the next event.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.nexteventandtime-Tuple{Dict}","page":"Public API","title":"DenseGillespieAlgorithm.nexteventandtime","text":"nexteventandtime(rates::Dict)\n\nSample a exponential distributed random variable to determine the time for the next event and calls `choose_event`.\nThe return value is a triple consiting of the envent index and trait returned by `choose_event` and the time to the next event.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.onestep!-Union{Tuple{F2}, Tuple{F1}, Tuple{Any, Any, Any, Any, Any, F1, F2}} where {F1, F2}","page":"Public API","title":"DenseGillespieAlgorithm.onestep!","text":"onestep!(x_0,rates,t_0,t_end,par,ex!::F1,r!::F2)\n\nExecute one step of the evolution by modifying `x_0` and `rates` and returning the current time `t_0`.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.run_gillespie!-Union{Tuple{F3}, Tuple{F2}, Tuple{F1}, Tuple{Any, Any, Any, F1, F2, Any, Any}} where {F1, F2, F3}","page":"Public API","title":"DenseGillespieAlgorithm.run_gillespie!","text":"    run_gillespie!(time,n₀,par,execute!,rates!,initrates,population_history[,hstart=0,statistic!])\n\nRun a exact stochastic simulation, return and fill the population_history.\n\nArguments\n\ntime::AbstracVector: time interval for the simulation\nn₀: initial population state\npar: additional parameter (gets passed to execute! and rates!)\nexecute!: execute function\nrates!: rates function\ninitrates: initial rates\npopulation_history: empty population history\nhstart=0: time shift for parameter change (opitonal)\nstatistic!: additional statistic function (optional)\n\nExtended help\n\nNote that n₀,initrates,population_history all three get modified during the simulation\nThe algorithm expects the execute! function to have the following signature   julia   execute!(i::Number,n₀,par)   where the i is the event that gets executed and the population state n₀ gets modified accordingly.   The only exception is when the initrates are given as a dictionary. In that case the signature is execute!(i,trait,n₀,initrates,par), where 'trait'   is the key that is modified.\nThe algorithm expects the rates! function to have the following signature   julia   rates!(initrates,n₀,par)   where the rates get modified according to the current population state given in n₀.\nThe algorithm expects the statistic! function to have the following signature   julia   statistic!(population_hist,t,n₀,par)   where the population history gets modified at position t with the current population state n₀.\nNote that the population_history needs to be accessable via index from 1 to length(time), or if hstart is given from 1+hstart to length(time)+hstart. Unless a specified statistic! function is given.\nNote that the initial population state n₀ must match the population_history in the sense that population_history :: Vector{typeof(n₀)}.  Unless a specified statistic! function is given.\nThe parameter variable par is passed through all functions (execute!,rates!,statistics!), thereby affording the user additional flexibility.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.saveonestep!-Tuple{Any, Any, Dict{<:Any, <:Number}, Any}","page":"Public API","title":"DenseGillespieAlgorithm.saveonestep!","text":" saveonestep!(pop_hist,index,ps,par)\n\n Sav one step of the simulation. Generic method if no explicit statistic! function is given.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.stop!-Union{Tuple{F1}, Tuple{Any, Any, Any, Any, F1}} where F1","page":"Public API","title":"DenseGillespieAlgorithm.stop!","text":"stop!(pop_hist,index,n0,par,stat!)\n\nFill the remaining population history with the (statistic of) the current population state if the evolution came to a halt.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.sumsumdict-Tuple{Any}","page":"Public API","title":"DenseGillespieAlgorithm.sumsumdict","text":"sumsumdict(D::Dict{String,Vector})\n\nCalculate the sum of the sums of the vectors that are the values of a dictionary.\n\n\n\n\n\n","category":"method"},{"location":"DGAPackage.html#The-DenseGillespieAlgorithm","page":"Home","title":"The DenseGillespieAlgorithm","text":"","category":"section"},{"location":"DGAPackage.html","page":"Home","title":"Home","text":"This package implements a version of the Gillespies algorithm that performs exact stochastic simulations for dense problems. The Gillespie algorithm[Gillespie76], introduced by Daniel Gillespie in 1976, is a fundamental tool for simulating the time evolution of systems with discrete, stochastic events, particularly in contexts like biochemical reactions and population dynamics. Its applications are particularly prevalent in contexts such as biochemical reactions and population dynamics. The Gillespie Algorithm is employed to simulate the behaviour of systems wherein reactions or events occur at random intervals. The algorithm generates a sequence of events and their timings by first calculating the rates at which different events or reactions occur. Subsequently, the time until the next event is determined based on these rates, and the type of event that occurs next is selected according to its probability. In the final step, the system state is updated based on the event, and the process is repeated.","category":"page"},{"location":"DGAPackage.html","page":"Home","title":"Home","text":"[Gillespie76]: D.T. Gillespie. A general method for numerically simulating the stochastic time evolution of coupled chemical reactions. Journal of Computational Physics, 22(4):403-434, 1976","category":"page"},{"location":"DGAPackage.html","page":"Home","title":"Home","text":"The Gillespie algorithm is a highly renowned and widely utilised technique across diverse communities and ecosystems. A particularly efficient, flexible and comprehensive implementation can be found in the JumpProcess.jl package within the SciML ecosystem. We strongly recommend the use of this framework wherever feasible. ","category":"page"},{"location":"DGAPackage.html","page":"Home","title":"Home","text":"However, the majority of implementations of the Gillespie algorithm require prior knowledge of all potential types and all reactions between these types before the reaction commences. A classic illustration of this is the SIR model (see Examples. The objective of our implementation in this package is to eliminate this restriction and permit the consideration of both high-dimensional systems, where the precise interactions between every conceivable combination are theoretically possible but practically infeasible, and additionally, systems where the trait space is uncountable, such as the real line. In both cases, the number of distinct traits that are present at any given time is finite, given that the population size is limited. However, new types emerge during the course of the simulation, and the interactions between these types are determined by their specific characteristics. ","category":"page"},{"location":"DGAPackage.html#Manual-Outline","page":"Home","title":"Manual Outline","text":"","category":"section"},{"location":"DGAPackage.html","page":"Home","title":"Home","text":"Pages = [\"manual.md\",\"examples.md\",\"perform.md\",\"index.md\"]","category":"page"},{"location":"DGAPackage.html#Index","page":"Home","title":"Index","text":"","category":"section"},{"location":"DGAPackage.html","page":"Home","title":"Home","text":"Modules = [DenseGillespieAlgorithm]\nOrder   = [:constant, :type, :function, :macro]","category":"page"}]
}