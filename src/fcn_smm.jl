@everywhere using SharedArrays
@everywhere include("fcn_moments.jl")
@everywhere include("fcn_results.jl")
include("list_models.jl")

"""
    smm(data, setup, mom_set, seed)

For a matrix of observations, perform parallel SMM estimation.

# Arguments
- `data::Array{Float64,2}`: matrix of observations
- `setup::Dict`: full setup dictionary
- `mom_set::Array`: moment set
- `seed::Int`: seed number
"""
function smm(data::Array{Float64,2}, setup::Dict, mom_set::Array, seed::Int)
    n_par = length(get_model_cali(setup["mod"])) # number of parameters
    n_rep = setup["smm"]["rep"] # number of repetitions

    # initialize shared arrays
    results_par = SharedArray{Float64}(n_par, n_rep) # estimated parameters
    results_j = SharedArray{Float64}(n_rep) # J-values

    # distribute work among workers
    @sync @distributed for i in 1:n_rep
        Random.seed!(10000*i+seed) # set random number generator seed

        # estimate parameters for the i-th column of the matrix of observations
        results_par[:,i], results_j[i] = gen_results(data[:,i], setup, mom_set)
    end

    return (Array(results_par), mom_set, Array(results_j))
end
