using Statistics

"""
    get_rmse(results, cali)

Calculate RMSE values for an array of results sets.

# Arguments
- `results::Array`: array of results sets
- `cali::Array{Float64}`: model calibration
"""
function get_rmse(results::Array, cali::Array{Float64})
    n_par = length(cali) # number of parameters
    n_res = length(results) # number of results sets

    rmse = zeros(n_par, n_res) # matrix of RMSE values

    for i in 1:n_res
        for j in 1:n_par
            # calculate bias and variance for a parameter in a results set
            bias_cur = abs(mean(results[i][1][j,:]) - cali[j])
            var_cur = std(results[i][1][j,:])^2

            # calculate normalized RMSE value
            rmse[j,i] = sqrt(bias_cur^2 + var_cur) / abs(cali[j])
        end
    end

    # take mean of normalized RMSE values for each results set
    rmse_vec = [mean(rmse[:,i]) for i in 1:n_res]

    return rmse_vec
end
