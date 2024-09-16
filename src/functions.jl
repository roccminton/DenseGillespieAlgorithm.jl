"""
        run_gillespie(time::AbstractVector,n₀,par::Parameter)

Run a exact stochastic simulation over the Interval time with an initial Population
`n₀`, Parameter Set `par` and Model Configuration `conf`. The output of the function
is a Vector{PopulationState} of length `length(time)` with the history of the
states of the population during the simulation.
"""

function run_gillespie!(time,n₀,par,execute!::F1,rates!::F2,initrates,population_history;hstart=0,statistic!::F3 = saveonestep!) where {F1,F2,F3}

    mainiteration!(
        population_history,
        initrates,
        n₀,
        convert(Float64,time[1]),
        time,
        par,
        execute!,
        rates!,
        statistic!,
        hstart
    )

end

"""
    mainiteration!

    For OneType model where the population state is a number and
    the population history is a vector.
"""
function mainiteration!(pop_hist,rates,n0::Real,ct,time,par,ex!::F1,r!::F2,stat!::F3,hstart) where {F1,F2,F3}
    #run simulation
    @showprogress for (index,step) in enumerate(time)
        #move index
        index += hstart
        #save one step evolution
        stat!(pop_hist,index,n0,par)
        #execute one step of the simulation
        ct, n0 = onestep!(n0,rates,ct,step,par,ex!,r!)
        #check if step was completed or evolution stopped inbetween
        @fastmath step-ct > 0.0 && break
    end
end

"""
    mainiteration!

    For models where population states are dictionaries (or vectors) with the traits as
    keys and the subpopulation size as values. Here the population history is
    itselfe a dictionary with the same keys and the individual subpopulation
    history as a vector for every trait.
"""
function mainiteration!(pop_hist,rates,n0,ct,time,par,ex!::F1,r!::F2,stat!::F3,hstart) where {F1,F2,F3}
    #run simulation
    @showprogress for (index,step) in enumerate(time)
        #move index
        index += hstart
        #clean up dictionary by deleting keys with value zero
        #dropzeros!(n0)
        #save one step evolution
        stat!(pop_hist,index,n0,par)
        #execute one step of the simulation
        ct = onestep!(n0,rates,ct,step,par,ex!,r!)
        #check if step was completed or evolution stopped inbetween
        @fastmath step-ct > 0.0 && (stop!(pop_hist,index,n0,par,stat!); break)
    end
end

historylength(population_history::Vector,par) = length(population_history)
historylength(population_history::Matrix,par) = length(view(population_history,:,1))
historylength(population_history,par) = par.historylength

function dropzeros!(ps::Dict{<:Any,<:Number})
    for (x,nₓ) ∈ ps
        iszero(nₓ) && delete!(ps,x)
    end
    nothing
end

function dropzeros!(ps::Dict{<:Any,<:Vector})
    for (x,vₓ) ∈ ps
        iszero(vₓ[1]) && delete!(ps,x)
    end
    nothing
end

dropzeros!(ps) = nothing

function saveonestep!(pop_hist,index,ps::Dict{<:Any,<:Number},par)
    for (x,nₓ) in ps
        #if the trait x is new to the population setup an empty history
        !haskey(pop_hist,x) && (pop_hist[x] = spzeros(valtype(ps),historylength(pop_hist,par)))
        #add the population history for the given trait
        pop_hist[x][index] = nₓ
    end
end

function saveonestep!(pop_hist,index,ps::Dict{<:Any,<:Vector},par)
    for (x,vₓ) in ps
        #if the trait x is new to the population setup an empty history
        !haskey(pop_hist,x) && (pop_hist[x] = zeros(eltype(valtype(pop_hist)),historylength(pop_hist,par)))
        #add the population history for the given trait
        pop_hist[x][index] = vₓ[1]
    end
end

function saveonestep!(pop_hist,index,ps,par)
    view(pop_hist,index,:) .= ps
end

function stop!(pop_hist,index,n0,par,stat!::F1) where {F1}
    for i ∈ index+1:historylength(pop_hist,par)
        stat!(pop_hist,i,n0,par)
    end
    nothing
end

"""
Executes one step of the evolution by modifying `x_0` and `rates`.
"""
function onestep!(x_0,rates,t_0,t_end,par,ex!::F1,r!::F2) where {F1,F2}
    while t_0 ≤ t_end
        r!(rates,x_0,par)
        #choose next event and event time
        i, dt = nexteventandtime(rates)
        #Population in absorbing state
        iszero(i) && break
        #update time
        t_0 += dt
        #execute event
        ex!(i,x_0,par)
    end
    return t_0
end

function onestep!(x_0,rates::Dict,t_0,t_end,par,ex!::F1,r!::F2) where {F1,F2}
    while t_0 ≤ t_end
        r!(rates,x_0,par)
        #choose next event and event time
        i, trait, dt = nexteventandtime(rates)
        #Population in absorbing state
        iszero(i) && break
        #update time
        t_0 += dt
        #execute event
        ex!(i,trait,x_0,rates,par)
    end
    return t_0
end

function onestep!(x_0::Real,rates::Vector,t_0,t_end,par,ex!::F1,r!::F2) where {F1,F2}
    while t_0 ≤ t_end
        r!(rates,x_0,par)
        #choose next event and event time
        i, dt = nexteventandtime(rates)
        #Population in absorbing state
        iszero(i) && break
        #update time
        t_0 += dt
        #execute event
        x_0 = ex!(i,x_0,par)
    end
    return t_0, x_0
end

function onestep!(x_0::Dict,rates::Vector,t_0,t_end,par,ex!::F1,r!::F2) where {F1,F2}
    while t_0 ≤ t_end
        r!(rates,x_0,par)
        #choose next event and event time
        i, dt = nexteventandtime(rates)
        #Population in absorbing state
        iszero(i) && break
        #update time
        t_0 += dt
        #execute event
        ex!(i,x_0,par)
    end
    return t_0
end

"""
    nexteventandtime(rates::Vector{Float64})

Samples a exponential distributed random variable to determine the time for the next
event and calls `choose_event`. The return value is a tuple consiting of the
envent index returned by `choose_event` and the time to the next event.
"""
function nexteventandtime(rates)
    #calculate total event rate
    @fastmath total_rate = sum(rates)
    #Population in absorbing state if sum of rates is zero
    iszero(total_rate) && return 0, 0.0
    #sample next event time
    dt = -(1/total_rate)*log(rand())
    #choose event
    i = chooseevent(rates,total_rate)
    return i, dt
end

function nexteventandtime(rates::Dict)
    #calculate total event rate
    total_rate = sumsumdict(rates)
    #Population in absorbing state if sum of rates is zero
    iszero(total_rate) && return 0, zero(keytype(rates)), 0.0
    #sample next event time
    dt = -(1/total_rate)*log(rand())
    #choose event
    i, trait = chooseevent(rates,total_rate)
    return i, trait, dt
end

"""
    chooseevent(rates::Vector{Float64},total_rate::Float64;<keyword arguments>)

From the vector of total rates choose at random one of the indices of the vector
according to their rates.
The value 0 is returned if the total rates are positive, but too smale to let
the evolution continue. The maximum number of tries is set by `max_try=1000`.
"""
function chooseevent(rates,total_rate)
    #make it a uniform random variable in (0,total_rate)
    rndm = rand()*total_rate
    #choose the rate at random
    @inbounds for (index, rate) in enumerate(rates)
        rndm -= rate
        rndm ≤ 0.0 && return index
    end
    return 0
end

function chooseevent(rates::Dict,total_rate)
    #make it a uniform random variable in (0,total_rate)
    rndm = rand()*total_rate
    #choose the rate at random
    for (trait, rate) in rates
        for (index,r) in enumerate(rate)
            rndm -= r
            rndm ≤ 0.0 && return index, trait
        end
    end
    return 0, zero(keytype(rates))
end

"""
For a dictionary with vectors of values calculates the sum of all
values of all vectors combined.
"""
function sumsumdict(D)
    result = zero(eltype(valtype(D)))
    for v in values(D)
        result += sum(v)
    end
    return result
end
