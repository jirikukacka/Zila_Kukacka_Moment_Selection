using JLD

include("list_models.jl")
include("fcn_watson.jl")

"""
    gen_series(model, obs, burn, theta)

Generate an empirical series using the specified model and parameters.

# Arguments
- `mod_set::Dict`: model setup dictionary
- `obs::Int`: number of observations
- `burn::Int`: burn-in period length
- `theta::Array{Float64}`: parameter values
"""
function gen_series(mod_set::Dict, obs::Int, burn::Int, theta::Array{Float64})
    # retrieve model's implementation and generate data
    model = get_model_func(mod_set)
    data = model(obs, burn, theta)

    return data
end


"""
    gen_data(setup)

Generate a dataset based on the setup.

# Arguments
- `setup::Dict`: full setup dictionary
"""
function gen_data(setup::Dict)
    smm_rep = setup["smm"]["rep"] # number of repetition
    mod_obs = setup["mod"]["obs"] # number of observations
    mod_burn = setup["mod"]["burn"] # burn-in period length

    # initialize array to store data
    data = zeros(mod_obs, smm_rep)

    mod_cali = get_model_cali(setup["mod"]) # model's calibration

    # produce a pseudo-empirical dataset
    for i in 1:smm_rep
        data[:,i] = gen_series(setup["mod"], mod_obs, mod_burn, mod_cali)
    end

    return data
end


"""
    load_data(setup)

Load data from the file specified in the setup.

# Arguments
- `setup::Dict`: full setup dictionary
"""
function load_data(setup::Dict)
    # load vector of empirical observations from JLD file and duplicate it
    emp_series = JLD.load(datadir(setup["smm"]["emp"]), "data")
    data = repeat(emp_series, 1, setup["smm"]["rep"])

    return data
end