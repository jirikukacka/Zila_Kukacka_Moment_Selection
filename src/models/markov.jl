using Distributions
using Random

"""
    markov_chain(obs, m0, gamma)

Generate one realization of length `obs` of a Markov chain with two states valued `m0` and `2-m0`.

# Arguments
- `obs::Int`: number of observations
- `m0::Float64`: multiplier value
- `gamma::Float64`: renewal probability
"""
function markov_chain(obs::Int, m0::Float64, gamma::Float64)
    ret = zeros(obs) # array of realized values
    ret_last = rand() > 0.5 ? m0 : 2-m0 # last realized value

    for i in 1:obs
        # state is renewed with probability gamma
        if rand() < gamma
            # if state was renewed, choose the new state with equal probability
            if rand() > 0.5
                ret_last = 2-ret_last # update realized value if state has changed
            end
        end
        ret[i] = ret_last # save realized value
    end

    return ret
end


"""
    markov(obs, burn, theta, k)

Generate data from Markov-switching model with 2^k states.

# Arguments
- `obs::Int`: number of observations
- `burn::Int`: burn-in period length
- `theta::Array{Float64}`: parameter values
- `k::Int`: number of states factor
"""
function markov(obs::Int, burn::Int, theta::Array{Float64}, k::Int)
    # estimated parameters
    m0 = theta[1] # multiplier value
    sigma = theta[2] # constant scale factor

    # simple calculations
    total_obs = obs+burn

    # data structures
    M = zeros(total_obs, k) # array of multipliers

    # generate shocks
    u = rand(Normal(0, 1), total_obs)

    # switching mechanism
    gamma = [2.0^-(k-i) for i in 1:k] # renewal probabilities

    # simulate individual Markov chains
    M = [markov_chain(total_obs, m0, gamma[i]) for i in 1:k]
    M = reduce(hcat, M)

    # calculate returns
    ret = sigma * sqrt.(prod(M, dims=2)) .* u # r_t = σ (∏_{i=1}^k M_t^{(i)}) u_t

    # discard burn-in period
    ret = ret[(burn+1):total_obs]

    return ret
end


"""
markov4state(obs, burn, theta)

Generate data from Markov-switching model with 4 states.

# Arguments
- `obs::Int`: number of observations
- `burn::Int`: burn-in period length
- `theta::Array{Float64}`: parameter values
"""
function markov4state(obs::Int, burn::Int, theta::Array{Float64})
    return markov(obs, burn, theta, 2)
end


"""
markov16state(obs, burn, theta)

Generate data from Markov-switching model with 16 states.

# Arguments
- `obs::Int`: number of observations
- `burn::Int`: burn-in period length
- `theta::Array{Float64}`: parameter values
"""
function markov16state(obs::Int, burn::Int, theta::Array{Float64})
    return markov(obs, burn, theta, 4)
end


"""
markov256state(obs, burn, theta)

Generate data from Markov-switching model with 256 states.

# Arguments
- `obs::Int`: number of observations
- `burn::Int`: burn-in period length
- `theta::Array{Float64}`: parameter values
"""
function markov256state(obs::Int, burn::Int, theta::Array{Float64})
    return markov(obs, burn, theta, 8)
end
