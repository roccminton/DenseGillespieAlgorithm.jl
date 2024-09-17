var documenterSearchIndex = {"docs":
[{"location":"index.html#The-DenseGillespieAlgorithm-Module","page":"Index","title":"The DenseGillespieAlgorithm Module","text":"","category":"section"},{"location":"index.html#Module-Index","page":"Index","title":"Module Index","text":"","category":"section"},{"location":"index.html","page":"Index","title":"Index","text":"Modules = [DenseGillespieAlgorithm]\nOrder   = [:constant, :type, :function, :macro]","category":"page"},{"location":"index.html#Detailed-API","page":"Index","title":"Detailed API","text":"","category":"section"},{"location":"index.html","page":"Index","title":"Index","text":"Modules = [DenseGillespieAlgorithm]\nOrder   = [:constant, :type, :function, :macro]","category":"page"},{"location":"index.html#DenseGillespieAlgorithm.chooseevent-Tuple{Any, Any}","page":"Index","title":"DenseGillespieAlgorithm.chooseevent","text":"chooseevent(rates::Vector{Float64},total_rate::Float64)\n\nChoose from the vector of total rates at random one of the indices of the vector according to their rates.\nThe value 0 is returned if the total rates are positive, but too smale to let the evolution continue.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.chooseevent-Tuple{Dict, Any}","page":"Index","title":"DenseGillespieAlgorithm.chooseevent","text":"chooseevent(rates::Dict,total_rate::Float64)\n\nChoose from the dictionary of total rates at random one of the keys of the dictionary according to their values.\nThe value 0 is returned if the total rates are positive, but too smale to let the evolution continue.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.dropzeros!-Tuple{Any}","page":"Index","title":"DenseGillespieAlgorithm.dropzeros!","text":"dropzeros!(ps)\n\nDo nothing for non-dictionary inputs.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.dropzeros!-Tuple{Dict{<:Any, <:Vector}}","page":"Index","title":"DenseGillespieAlgorithm.dropzeros!","text":"dropzeros!(ps::Dict{Any,Vector})\n\nEliminates all key value pairs for which the the firts entry of the vector of the value is zero.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.historylength-Tuple{Vector, Any}","page":"Index","title":"DenseGillespieAlgorithm.historylength","text":"historylength(population_hisotry,par)\n\nReturn the simulation time based on the length of the population history.\nIf the population history is neither a `Vector` nor a `Matrix` it is assumed that the Parameter has a field called `historylength` that is then returned.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.mainiteration!-Union{Tuple{F3}, Tuple{F2}, Tuple{F1}, Tuple{Any, Any, Any, Any, Any, Any, F1, F2, F3, Any}} where {F1, F2, F3}","page":"Index","title":"DenseGillespieAlgorithm.mainiteration!","text":"mainiteration!(pop_hist,rates,n0,ct,time,par,ex!::F1,r!::F2,stat!::F3,hstart)\n\nMainiteration of the GillespieAlgorithm for complex models.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.mainiteration!-Union{Tuple{F3}, Tuple{F2}, Tuple{F1}, Tuple{Any, Any, Real, Any, Any, Any, F1, F2, F3, Any}} where {F1, F2, F3}","page":"Index","title":"DenseGillespieAlgorithm.mainiteration!","text":"mainiteration!(pop_hist,rates,n0::Real,ct,time,par,ex!::F1,r!::F2,stat!::F3,hstart)\n\nMainiteration of the GillespieAlgorithm for OneType Models where the population state is a number and the population history is a vector.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.nexteventandtime-Tuple{Any}","page":"Index","title":"DenseGillespieAlgorithm.nexteventandtime","text":"nexteventandtime(rates::Vector{Float64})\n\nSample a exponential distributed random variable to determine the time for the next event and calls `choose_event`.\nThe return value is a tuple consiting of the envent index returned by `choose_event` and the time to the next event.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.nexteventandtime-Tuple{Dict}","page":"Index","title":"DenseGillespieAlgorithm.nexteventandtime","text":"nexteventandtime(rates::Dict)\n\nSample a exponential distributed random variable to determine the time for the next event and calls `choose_event`.\nThe return value is a triple consiting of the envent index and trait returned by `choose_event` and the time to the next event.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.onestep!-Union{Tuple{F2}, Tuple{F1}, Tuple{Any, Any, Any, Any, Any, F1, F2}} where {F1, F2}","page":"Index","title":"DenseGillespieAlgorithm.onestep!","text":"onestep!(x_0,rates,t_0,t_end,par,ex!::F1,r!::F2)\n\nExecute one step of the evolution by modifying `x_0` and `rates` and returning the current time `t_0`.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.run_gillespie!-Union{Tuple{F3}, Tuple{F2}, Tuple{F1}, Tuple{Any, Any, Any, F1, F2, Any, Any}} where {F1, F2, F3}","page":"Index","title":"DenseGillespieAlgorithm.run_gillespie!","text":"    run_gillespie!(time,n₀,par,execute!,rates!,initrates,population_history[,hstart=0,statistic!])\n\nRun a exact stochastic simulation, return and fill the population_history.\n\n#Arguments\n\ntime::AbstracVector: time interval for the simulation\nn₀: initial population state\npar: additional parameter (gets passed to execute! and rates!)\nexecute!: execute function\nrates!: rates function\ninitrates: initial rates\npopulation_history: empty population history\nhstart=0: time shift for parameter change (opitonal)\nstatistic!: additional statistic function (optional)\n\nExtended help\n\nNote that n₀,initrates,population_history all three get modified during the simulation\nThe algorithm expects the execute! function to have the following signature   julia   execute!(i::Number,n₀,par)   where the i is the event that gets executed and the population state n₀ gets modified accordingly.   The only exception is when the initrates are given as a dictionary. In that case the signature is execute!(i,trait,n₀,initrates,par), where 'trait'   is the key that is modified.\nThe algorithm expects the rates! function to have the following signature   julia   rates!(initrates,n₀,par)   where the rates get modified according to the current population state given in n₀.\nThe algorithm expects the statistic! function to have the following signature   julia   statistic!(population_hist,t,n₀,par)   where the population history gets modified at position t with the current population state n₀.\nNote that the population_history needs to be accessable via index from 1 to length(time), or if hstart is given from 1+hstart to length(time)+hstart. Unless a specified statistic! function is given.\nNote that the initial population state n₀ must match the population_history in the sense that population_history :: Vector{typeof(n₀)}.  Unless a specified statistic! function is given.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.saveonestep!-Tuple{Any, Any, Dict{<:Any, <:Number}, Any}","page":"Index","title":"DenseGillespieAlgorithm.saveonestep!","text":" saveonestep!(pop_hist,index,ps,par)\n\n Sav one step of the simulation. Generic method if no explicit statistic! function is given.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.stop!-Union{Tuple{F1}, Tuple{Any, Any, Any, Any, F1}} where F1","page":"Index","title":"DenseGillespieAlgorithm.stop!","text":"stop!(pop_hist,index,n0,par,stat!)\n\nFill the remaining population history with the (statistic of) the current population state if the evolution came to a halt.\n\n\n\n\n\n","category":"method"},{"location":"index.html#DenseGillespieAlgorithm.sumsumdict-Tuple{Any}","page":"Index","title":"DenseGillespieAlgorithm.sumsumdict","text":"sumsumdict(D::Dict{String,Vector})\n\nCalculate the sum of the sums of the vectors that are the values of a dictionary.\n\n\n\n\n\n","category":"method"}]
}
