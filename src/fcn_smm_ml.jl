include("fcn_smm.jl")
include("fcn_data.jl")
include("fcn_metric.jl")
include("fcn_utils.jl")
include("fcn_watson.jl")
include("list_models.jl")

"""
    smm_init(setup)

Prepare dataset and initialize selected estimation procedure.

If `setup["smm"]["emp"]` is not set, generate pseudo-empirical data. Otherwise,
load empirical data from the specified file located in the `data` folder.

# Arguments
- `setup::Dict`: full setup dictionary
"""
function smm_init(setup::Dict)
    Random.seed!(42) # set random seed

    if setup["smm"]["emp"] isa Nothing
        data = gen_data(setup) # generate pseudo-empirical data
    else
        data = load_data(setup) # load empirical data
        setup["mod"]["obs"] = size(data, 1) # adjust length of generated series
    end

    # generate name for results folder
    foldername = make_foldername(setup["ml"], setup["mod"], setup["opt"], setup["wgt"])

    # forward stepwise moment selection
    if setup["ml"]["method"] == "fsw"
        results = fsw(data, foldername, setup)
    # backward stepwise moment elimination
    elseif setup["ml"]["method"] == "bsw"
        results = bsw(data, foldername, setup)
    # single moment set estimation
    elseif setup["ml"]["method"] == "sms"
        results = sms(data, foldername, setup)
    else
        error("Underfined method of estimation.")
    end

    return results
end


"""
    fsw(data, setup)

Apply forward stepwise moment selection procedure to SMM.

# Arguments
- `data::Array{Float64,2}`: matrix of observations
- `foldername::String`: name for results folder
- `setup::Dict`: full setup dictionary
"""
function fsw(data::Array{Float64,2}, foldername::String, setup::Dict)
    results = [] # array of estimation results

    model_cali = get_model_cali(setup["mod"]) # model calibration

    mom_set = repeat([false], length(SET["full"])) # initialize moment set

    # forward stepwise moment selection
    while sum(mom_set) < length(mom_set)
        flushln("Round $(sum(mom_set)+1) of forward stepwise moment selection.")

        results_cur = [] # array of estimation results for the current round
        
        # attempt to add each unused moment to the moment set
        for i in findall(iszero, mom_set)
            # add moment i to the moment set
            mom_set_cur = copy(mom_set)
            mom_set_cur[i] = true

            # estimate parameters for the expanded moment set
            flushln("Producing results for:\n$mom_set_cur.")
            res = smm(data, setup, mom_set_cur, length(results))

            push!(results_cur, res)
            push!(results, res)
        end

        # save estimation results to disk
        save_results(results_cur, foldername, length(results)-length(results_cur))

        # calculate RMSE for each moment set expansion
        results_metric = get_rmse(results_cur, model_cali)

        # choose moment set expansion with lowest RMSE
        mom_set = copy(results_cur[argmin(results_metric)][2])
        flushln("The best selection of moments was chosen to be:\n$mom_set.")
    end

    return results
end


"""
    bsw(data, setup)

Apply backward stepwise moment elimination procedure to SMM.

# Arguments
- `data::Array{Float64,2}`: matrix of observations
- `foldername::String`: name for results folder
- `setup::Dict`: full setup dictionary
"""
function bsw(data::Array{Float64,2}, foldername::String, setup::Dict)
    results = [] # array of estimation results

    model_cali = get_model_cali(setup["mod"]) # model calibration

    mom_set = repeat([true], length(SET["full"])) # initialize moment set

    # backward stepwise moment elimination
    while sum(mom_set) > 1
        flushln("Round $(length(mom_set)-sum(mom_set)+1) of backward stepwise moment elimination.")

        results_cur = [] # array of estimation results for the current round

        # attempt to remove each used moment from the moment set
        for i in findall(isone, mom_set)
            # remove moment i from the moment set
            mom_set_cur = copy(mom_set)
            mom_set_cur[i] = false

            # estimate parameters for the reduced moment set
            flushln("Producing results for:\n$mom_set_cur.")
            res = smm(data, setup, mom_set_cur, length(results))

            push!(results_cur, res)
            push!(results, res)
        end

        # save estimation results to disk
        save_results(results_cur, foldername, length(results)-length(results_cur))

        # calculate RMSE for each moment set reduction
        results_metric = get_rmse(results_cur, model_cali)

        # choose moment set reduction with lowest RMSE
        mom_set = copy(results_cur[argmin(results_metric)][2])
        flushln("The best selection of moments was chosen to be:\n$mom_set.")
    end

    return results
end


"""
    sms(data, setup)

Produce SMM estimation results for a single moment set.

# Arguments
- `data::Array{Float64,2}`: matrix of observations
- `foldername::String`: name for results folder
- `setup::Dict`: full setup dictionary
"""
function sms(data::Array{Float64,2}, foldername::String, setup::Dict)
    results = [] # array of estimation results

    if setup["ml"]["bench"] isa Nothing
        mom_set = setup["ml"]["set"]
    else
        mom_set = SET[setup["ml"]["bench"]]
    end

    # estimate parameters for the moment set
    flushln("Producing results for:\n$mom_set.")
    res = smm(data, setup, mom_set, length(results))
    
    push!(results, res)

    # save estimation results to disk
    save_results(results, foldername, 0)

    return results
end
