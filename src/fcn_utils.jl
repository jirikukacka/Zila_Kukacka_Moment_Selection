using JLD, Plots.PlotMeasures, StatsPlots

include("fcn_watson.jl")
include("fcn_metric.jl")
include("fcn_data.jl")
include("list_models.jl")

##############
# FOLDERNAME #
##############

"""
    make_foldername(setup)

Concatenate settings into a comprehensive foldername.

# Arguments
- `setup::Dict`: full setup dictionary
"""
function make_foldername(setup::Dict)
    alg = setup["ml"]["method"]
    if alg == "sms" && !isnothing(setup["ml"]["bench"])
        alg = setup["ml"]["bench"]
    end

    foldername = string(setup["mod"]["model"], "_",
                        alg, "_",
                        setup["mod"]["cali"], "_",
                        setup["mod"]["cons"], "_",
                        "rep", setup["smm"]["rep"], "_",
                        "obs", setup["mod"]["obs"], "_",
                        "burn", setup["mod"]["burn"], "_",
                        "inits", setup["opt"]["inits"], "_",
                        "sim", setup["opt"]["sim"], "_",
                        "iter", setup["opt"]["iter"], "_",
                        setup["wgt"]["method"], "_",
                        "blocksize", setup["wgt"]["blocksize"], "_",
                        "bootsize", setup["wgt"]["bootsize"], "_",
                        "blockcount", setup["wgt"]["blockcount"])

    return foldername
end

"""
    foldername_alg_replace(foldername, newalg)

Replace algorithm name in a foldername for a given string.

# Arguments
- `foldername::String`: original foldername
- `newalg::String`: updated foldername with a new algorithm string
"""
function foldername_alg_replace(foldername::String, newalg::String)
    foldername_split = split(foldername, "_")
    foldername_split[2] = newalg

    return join(foldername_split, "_")
end

"""
    foldername_mod_set_retrieve(foldername)

Retrieve model setup dictionary from from a given foldername.

# Arguments
- `foldername::String`: foldername to retrieve model setup from
"""
function foldername_mod_set_retrieve(foldername)
    foldername_split = split(foldername, "_")

    mod_set = Dict("model" => foldername_split[1],
                   "cali" => foldername_split[3],
                   "cons" => foldername_split[4])

    return mod_set
end


###########
# RESULTS #
###########

"""
    save_results(results, foldername, idx=0)

Save results at a given location.

# Arguments
- `results`: set of results
- `foldername::String`: location to store the set of results at
- `idx::Int`: initial indexing value
"""
function save_results(results, foldername::String, idx=0)
    for i in eachindex(results)
        save(resultsdir(foldername, "set$(idx+i).jld"), "results", results[i])
    end
end

"""
    load_results(foldername)

Load results from a given location.

# Arguments
- `foldername::String`: location to load the set of results from
"""
function load_results(foldername::String)
    folder_sets = [SubString(x, 1, 3) for x in readdir(resultsdir(foldername))]
    set_count = sum([x == "set" for x in folder_sets])

    results = []

    for i in 1:set_count
        content = load(resultsdir(foldername, "set$(i).jld"))
        try
            push!(results, (content["data"], content["model"]))
        catch
            push!(results, (content["results"]))
        end
    end

    return results
end

"""
    load_bench(foldername)

Load benchmarks for a given set of results.

# Arguments
- `foldername::String`: location of the set of results
"""
function load_bench(foldername::String)
    results = []

    for bench in ["chl4", "chl15", "fw9"]
        try
            content = load(benchdir(foldername_alg_replace(foldername, bench), "set1.jld"))
            try
                push!(results, (content["data"], content["model"]))
            catch
                push!(results, (content["results"]))
            end
        catch
            println("Folder containing $bench benchmark is not available!")
        end
    end

    return results
end


#########
# OTHER #
#########

"""
    results_best(results, cali)

Retrieve the best performing moment set set for each round of the selection 
process given a set of results and its calibration.

# Arguments
- `results`: set of results
- `cali`: calibration of the given set of results
"""
function results_best(results, cali)
    rounds = [sum(results[i][2]) for i in eachindex(results)]

    best_idx = []
    for round in unique(rounds)
        idx_temp = findall(rounds .== round)
        results_temp = [results[i] for i in idx_temp]
        rmse_temp = get_rmse(results_temp, cali)

        push!(best_idx, idx_temp[argmin(rmse_temp)])
    end

    rsb = [results[i] for i in best_idx]

    return (rsb, best_idx)
end

"""
    print_evolution(results, type)

Print evolution of the selection process given a set of results and the 
direction of the selection process.

# Arguments
- `results`: set of results
- `type::String`: `"bsw"` for backward stepwise moemnt elimination, `"fsw"` for 
forward stepwise moment selection
"""
function print_evolution(results, type::String)
    cali = get_model_cali(mod_set)

    mom_names = ["RAW-VAR", "RAW-KURT", "RAW-AC1", "RAW-AC2", "RAW-AC3",
                 "ABS-MEAN", "ABS-HILL2.5", "ABS-HILL5", "ABS-AC1", "ABS-AC5", 
                 "ABS-AC10", "ABS-AC15", "ABS-AC20", "ABS-AC25", "ABS-AC50", 
                 "ABS-AC100", "SQR-AC1", "SQR-AC5", "SQR-AC10", "SQR-AC15", 
                 "SQR-AC20", "SQR-AC25"]

    rsb, _ = results_best(results, cali)
    mom_series = [] # storage of added/removed moments in each round

    if type == "fsw" || type == "bsw"
        is_fsw = type == "fsw"

        last_best = is_fsw ? zeros(length(mom_names)) : ones(length(mom_names))

        for i in eachindex(rsb)
            cur_mom = rsb[i][2]
            mom_idx = is_fsw ? argmax(cur_mom - last_best) : mom_idx = argmax(last_best - cur_mom)
            last_best = cur_mom

            push!(mom_series, mom_names[mom_idx])
        end

        if !is_fsw
            push!(mom_series, mom_names[argmax(rsb[length(rsb)][2])])
        end
        
        ret = join(mom_series, " \$\\rightarrow\$ ")

        print(ret)
    else
        error("Undefined type inserted! Defined types are fsw and bsw.")
    end
end
