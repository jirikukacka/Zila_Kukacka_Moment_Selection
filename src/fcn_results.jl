using BlackBoxOptim

include("fcn_data.jl")
include("fcn_moments.jl")
include("fcn_weights.jl")
include("list_models.jl")

"""
    J_fcn(theta, setup, moments_emp, mom_set, weights, simlen)

Calculate the J value computed as the weighted difference between empirical 
moments and average simulated moments.

# Arguments
- `theta::Array{Float64}`: parameter values
- `setup::Dict`: full setup dictionary
- `moments_emp::Array`: array of empirical moments
- `mom_set::Array`: moment set
- `weights`: SMM weighting matrix
- `simlen::Int`: length of the simulated series
"""
function J_fcn(theta::Array{Float64}, setup::Dict, moments_emp::Array, mom_set::Array, weights, simlen::Int)
    moments_sim = zeros(sum(mom_set), setup["opt"]["sim"]) # array of simulated moments

    for i in 1:setup["opt"]["sim"]
        # simulate time series and calculate simulated moments
        data = gen_series(setup["mod"], simlen, setup["mod"]["burn"], theta)
        moments_sim[:,i] = gen_moments_sel(data, mom_set)
    end

    # calculate J value
    moments_diff = mean(moments_sim, dims=2) - moments_emp
    obj = transpose(moments_diff) * weights * moments_diff 

    return obj[1]
end


"""
    gen_results(data, setup, mom_set)

Produce SMM estimation results for a vector of observations.

# Arguments
- `data::Array{Float64,1}`: array of observations
- `setup::Dict`: full setup dictionary
- `mom_set::Array`: moment set
"""
function gen_results(data::Array{Float64,1}, setup::Dict, mom_set::Array)
    iter = setup["opt"]["iter"] # number of iterations
    search_range = get_model_cons(setup["mod"]) # search constraints

    # data structures
    results_parm = zeros(length(search_range), setup["opt"]["inits"]) # matrix of estimated parameters
    results_j = zeros(setup["opt"]["inits"]) # array of J values

    weights = gen_weights(setup, data, mom_set) # generate SMM weighting matrix
    moments_emp = gen_moments_sel(data, mom_set) # calculate empirical moments 

    for i in 1:setup["opt"]["inits"]
        # select parameters such that J value is minimized
        optout = bboptimize(theta -> J_fcn(theta, setup, moments_emp, mom_set, weights, 
                                           length(data)*setup["smm"]["simfactor"]),
                            SearchRange = search_range,
                            Method = :adaptive_de_rand_1_bin_radiuslimited,
                            NumDimensions = length(search_range),
                            MaxFuncEvals = iter,
                            TraceMode = :silent)

        results_parm[:,i] = best_candidate(optout) # retrieve found parameters
        results_j[i] = best_fitness(optout) # retrieve corresponding J value
    end

    return (results_parm[:, argmin(results_j)], minimum(results_j))
end
