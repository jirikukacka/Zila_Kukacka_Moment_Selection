using JLD, Plots.PlotMeasures, StatsPlots

include("fcn_watson.jl")
include("fcn_metric.jl")
include("fcn_data.jl")
include("list_models.jl")


###############
# FOLDERNAMES #
###############

# concatenate settings into results foldername
function make_foldername(ml_set, mod_set, opt_set, wgt_set)
    alg = ml_set["method"]
    if alg == "sms" && !isnothing(ml_set["bench"])
        alg = ml_set["bench"]
    end

    foldername = string(mod_set["model"], "_",
                        alg, "_",
                        mod_set["cali"], "_",
                        mod_set["cons"], "_",
                        "rep", smm_set["rep"], "_",
                        "obs", mod_set["obs"], "_",
                        "burn", mod_set["burn"], "_",
                        "inits", opt_set["inits"], "_",
                        "sim", opt_set["sim"], "_",
                        "iter", opt_set["iter"], "_",
                        wgt_set["method"], "_",
                        "blocksize", wgt_set["blocksize"], "_",
                        "bootsize", wgt_set["bootsize"], "_",
                        "blockcount", wgt_set["blockcount"])

    return foldername
end


# retrieve model settings from foldername
function foldername_mod_set_retrieve(foldername)
    foldername_split = split(foldername, "_")

    mod_set = Dict("model" => foldername_split[1],
                   "cali" => foldername_split[3],
                   "cons" => foldername_split[4])

    return mod_set
end


# replace algorithm string in foldername for an alternative
function foldername_alg_replace(foldername, newalg)
    foldername_split = split(foldername, "_")
    foldername_split[2] = newalg

    return join(foldername_split, "_")
end


###########
# RESULTS #
###########

# save results to the "results/" folder
function save_results(results, foldername, idx=0)
    for i in eachindex(results)
        save(resultsdir(foldername, "set$(idx+i).jld"), "results", results[i])
    end
end


# load results from the "results/" folder
function load_results(foldername)
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


# load benchmarks from the "results/bench/" folder
function load_bench(foldername)
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


function print_rmse(results, mod_set, set=nothing)
    cali = get_model_cali(mod_set)
    results_rmse = get_rmse(results, cali)

    println("RMSE values:")
    if isnothing(set)
        for i in eachindex(results_rmse)
            println("set $i : $(results_rmse[i])")
        end
    else
        println("set $set : $(results_rmse[set])")
    end
end


function print_evolution(results, type)
    cali = get_model_cali(mod_set)

    mom_names = ["RAW-VAR", "RAW-KURT", "RAW-AC1", "RAW-AC2", "RAW-AC3",
                 "ABS-MEAN", "ABS-HILL2.5", "ABS-HILL5", "ABS-AC1", "ABS-AC5", 
                 "ABS-AC10", "ABS-AC15", "ABS-AC20", "ABS-AC25", "ABS-AC50", 
                 "ABS-AC100", "SQR-AC1", "SQR-AC5", "SQR-AC10", "SQR-AC15", 
                 "SQR-AC20", "SQR-AC25"]

    rsb, _ = results_best(results, cali)
    mom_series = [] # storage of added/removed moments

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
        
        #ret = is_fsw ? join(mom_series, " > ") : join(mom_series, " < ")
        ret = join(mom_series, " \$\\rightarrow\$ ")

        print(ret)
    else
        error("Undefined type inserted! Defined types are fsw and bsw.")
    end
end
