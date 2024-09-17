"""
        run_gillespie!(time,n₀,par,execute!,rates!,initrates,population_history[,hstart=0,statistic!])

Run a exact stochastic simulation, return and fill the `population_history`.

#Arguments
- `time::AbstracVector`: time interval for the simulation
- `n₀`: initial population state
- `par`: additional parameter (gets passed to `execute!` and `rates!`)
- `execute!`: execute function
- `rates!`: rates function
- `initrates`: initial rates
- `population_history`: empty population history
- `hstart=0`: time shift for parameter change _(opitonal)_
- `statistic!`: additional statistic function _(optional)_

# Extended help
- Note that `n₀,initrates,population_history` all three get modified during the simulation
- The algorithm expects the `execute!` function to have the following signature
    ```julia
    execute!(i::Number,n₀,par)
    ```
    where the `i` is the event that gets executed and the population state `n₀` gets modified accordingly.
    The only exception is when the `initrates` are given as a dictionary. In that case the signature is `execute!(i,trait,n₀,initrates,par)`, where 'trait'
    is the key that is modified.
- The algorithm expects the `rates!` function to have the following signature
    ```julia
    rates!(initrates,n₀,par)
    ```
    where the rates get modified according to the current population state given in `n₀`.
- The algorithm expects the `statistic!` function to have the following signature
    ```julia
    statistic!(population_hist,t,n₀,par)
    ```
    where the population history gets modified at position t with the current population state `n₀`.
- Note that the `population_history` needs to be accessable via index from 1 to `length(time)`, or if `hstart` is given from `1+hstart` to `length(time)+hstart`. Unless a specified `statistic!` function is given.
- Note that the initial population state `n₀` must match the `population_history` in the sense that `population_history :: Vector{typeof(n₀)}`.  Unless a specified `statistic!` function is given.
- The parameter variable `par` is passed through all functions (`execute!,rates!,statistics!`), thereby affording the user additional flexibility.
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
    mainiteration!(pop_hist,rates,n0,ct,time,par,ex!::F1,r!::F2,stat!::F3,hstart)

    Mainiteration of the GillespieAlgorithm for complex models.
"""
function mainiteration!(pop_hist,rates,n0,ct,time,par,ex!::F1,r!::F2,stat!::F3,hstart) where {F1,F2,F3}
    #run simulation
    ProgressMeter.@showprogress for (index,step) in enumerate(time)
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
"""
    mainiteration!(pop_hist,rates,n0::Real,ct,time,par,ex!::F1,r!::F2,stat!::F3,hstart)

    Mainiteration of the GillespieAlgorithm for OneType Models where the population state is a number and the population history is a vector.
"""
function mainiteration!(pop_hist,rates,n0::Real,ct,time,par,ex!::F1,r!::F2,stat!::F3,hstart) where {F1,F2,F3}
    #run simulation
    ProgressMeter.@showprogress for (index,step) in enumerate(time)
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
    historylength(population_hisotry,par)

    Return the simulation time based on the length of the population history.
    If the population history is neither a `Vector` nor a `Matrix` it is assumed that the Parameter has a field called `historylength` that is then returned.
"""
historylength(population_history::Vector,par) = length(population_history)
historylength(population_history::Matrix,par) = length(view(population_history,:,1))
historylength(population_history,par) = par.historylength

function dropzeros!(ps::Dict{<:Any,<:Number})
    for (x,nₓ) ∈ ps
        iszero(nₓ) && delete!(ps,x)
    end
    nothing
end

"""
    dropzeros!(ps::Dict{Any,Vector})

    Eliminates all key value pairs for which the the firts entry of the vector of the value is zero.
"""
function dropzeros!(ps::Dict{<:Any,<:Vector})
    for (x,vₓ) ∈ ps
        iszero(vₓ[1]) && delete!(ps,x)
    end
    nothing
end
"""
    dropzeros!(ps)

    Do nothing for non-dictionary inputs.
"""
dropzeros!(ps) = nothing

"""
     saveonestep!(pop_hist,index,ps,par)

     Sav one step of the simulation. Generic method if no explicit statistic! function is given.
"""
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

"""
    stop!(pop_hist,index,n0,par,stat!)

    Fill the remaining population history with the (statistic of) the current population state if the evolution came to a halt.
"""
function stop!(pop_hist,index,n0,par,stat!::F1) where {F1}
    for i ∈ index+1:historylength(pop_hist,par)
        stat!(pop_hist,i,n0,par)
    end
    nothing
end

"""
    onestep!(x_0,rates,t_0,t_end,par,ex!::F1,r!::F2)

    Execute one step of the evolution by modifying `x_0` and `rates` and returning the current time `t_0`.
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

    Sample a exponential distributed random variable to determine the time for the next event and calls `choose_event`.
    The return value is a tuple consiting of the envent index returned by `choose_event` and the time to the next event.
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

"""
    nexteventandtime(rates::Dict)

    Sample a exponential distributed random variable to determine the time for the next event and calls `choose_event`.
    The return value is a triple consiting of the envent index and trait returned by `choose_event` and the time to the next event.
"""
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
    chooseevent(rates::Vector{Float64},total_rate::Float64)

    Choose from the vector of total rates at random one of the indices of the vector according to their rates.
    The value 0 is returned if the total rates are positive, but too smale to let the evolution continue.
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

"""
    chooseevent(rates::Dict,total_rate::Float64)

    Choose from the dictionary of total rates at random one of the keys of the dictionary according to their values.
    The value 0 is returned if the total rates are positive, but too smale to let the evolution continue.
"""
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
    sumsumdict(D::Dict{String,Vector})

    Calculate the sum of the sums of the vectors that are the values of a dictionary.
"""
function sumsumdict(D)
    result = zero(eltype(valtype(D)))
    for v in values(D)
        result += sum(v)
    end
    return result
end
