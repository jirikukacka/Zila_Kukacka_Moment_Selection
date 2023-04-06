using JLD, Plots, Plots.PlotMeasures, StatsPlots

include("fcn_watson.jl")
include("fcn_metric.jl")
include("fcn_utils.jl")
include("fcn_data.jl")
include("list_models.jl")


##########################
# PLOTS FOR PRESENTATION #
##########################

function plot_histogram(results, mod_set, i, folder; name=nothing)
    res = results[i][1]
    filename = isnothing(name) ? "set$(i)_$(folder)_.pdf" : "$(name)_$(folder).pdf"

    cali = get_model_cali(mod_set)
    cons = get_model_cons(mod_set)
    gnam = get_model_gnam(mod_set)
    gdim = get_model_gdim(mod_set)
    gord = get_model_gord(mod_set)

    subplot_tot = gdim[1] * gdim[2]
    plots = Vector{Plots.Plot{Plots.GRBackend}}(undef, subplot_tot)

    for i = eachindex(cali)
        res_par = res[i,:]
        plots[i] = density(res_par,
                           legend = :none,
                           title = gnam[i],
                           color = :black,
                           framestyle = :border,
                           xlims = cons[i],
                           xtickfont = font(6),
                           ytickfont = font(6))
        plots[i] = vline!([cali[i]], color = :red)
        plots[i] = vline!([mean(res_par)], color = :black)
        plots[i] = vline!([quantile(res_par, 0.025),
                           quantile(res_par, 0.975)],
                          color = :black, line = :dash)
    end

    for i = (length(cali)+1):subplot_tot
        plots[i] = plot(border = :none)
    end

    graph_dims = (gdim[2]*200, gdim[1]*105)
    plot_res = get_resulting_plot(plots, gdim, graph_dims, gord)

    savefig(plot_res, plotsdir(filename))
end


function plot_rmse(results, mod_set, folder; bench=[], bench_names=[])
    filename = "rmse_$(folder).pdf"

    cali = get_model_cali(mod_set)

    results_full = deepcopy(results)
    for b in bench
        push!(results_full, b)
    end

    rmse = get_rmse(results_full, cali)
    rounds = [sum(results_full[i][2]) for i in eachindex(results_full)]

    annot = text.(["$i  " for i in 1:length(results)], 5, :right)
    mcol = [:grey for i in 1:length(results)]

    for b in bench_names
        push!(annot, text("  $b", 5, :left, :red))
        push!(mcol, :red)
    end

    plot = scatter(rounds, rmse,
                   series_annotations=annot,
                   markercolor=mcol,
                   legend=:none,
                   xlabel=L"number\ of\ moments",
                   ylabel=L"RMSE",
                   xticks=0:1:maximum(rounds),
                   size=(1000,750),
                   left_margin=4mm,
                   bottom_margin=5mm)
    savefig(plot, plotsdir(filename))
end


function plot_rmse_best(results, mod_set, folder; bench=[], bench_names=[])
    filename = "rmse_best_$(folder).pdf"

    cali = get_model_cali(mod_set)

    (rsb, best_idx) = result_best(results, cali)

    results_full = deepcopy(rsb)
    for b in bench
        push!(results_full, b)
    end

    rmse = get_rmse(results_full, cali)
    rounds = [sum(results_full[i][2]) for i in eachindex(results_full)]

    annot = text.(["$(best_idx[i])  " for i in eachindex(rsb)], 5, :right)
    mcol = [:grey for i in 1:length(rsb)]

    for b in bench_names
        push!(annot, text("  $b", 5, :left, :red))
        push!(mcol, :red)
    end

    plot = scatter(rounds, rmse,
                   series_annotations=annot,
                   markercolor=mcol,
                   legend=:none,
                   xlabel=L"number\ of\ moments",
                   ylabel=L"RMSE",
                   xticks=0:1:maximum(rounds),
                   size=(1000,750),
                   left_margin=4mm,
                   bottom_margin=5mm)
    savefig(plot, plotsdir(filename))
end


function plot_evolution_par(results, mod_set, folder; bench=[], bench_names=[], ylims=nothing)
    filename = "evolution_par_$(folder).pdf"

    cali = get_model_cali(mod_set)
    cons = get_model_cons(mod_set)
    gnam = get_model_gnam(mod_set)
    gdim = get_model_gdim(mod_set)
    gord = get_model_gord(mod_set)

    if !isnothing(ylims)
        cons = ylims
    end

    rsb, _ = results_best(results, cali)
    rounds = [sum(rsb[i][2]) for i in eachindex(rsb)]

    subplot_tot = gdim[1] * gdim[2]
    plots = Vector{Plots.Plot{Plots.GRBackend}}(undef, subplot_tot)

    for i = eachindex(cali)
        est = []
        ub = []
        lb = []

        for j in eachindex(rsb)
            est_cur = mean(rsb[j][1][i,:])
            push!(est, est_cur)

            b_cur = (1.960 * std(rsb[j][1][i,:])) #/ sqrt(size(rsb[1][1], 2))
            push!(ub, est_cur + b_cur)
            push!(lb, est_cur - b_cur)
        end

        plots[i] = plot(rounds, est,
                        legend=:none,
                        xticks=1:2:22,
                        title = gnam[i],
                        xlabel=L"number\ of\ moments",
                        xguidefontsize=9,
                        xlims=(0.5,22.5),
                        ylim=cons[i],
                        color=:black,
                        xtickfont = font(6),
                        ytickfont = font(6),
                        bottom_margin=2.5mm)
        
        #plot!(rounds, ub,
        #      color=:black,
        #      linestyle=:dash,
        #      linewidth=2)

        #plot!(rounds, lb,
        #      color=:black,
        #      linestyle=:dash,
        #      linewidth=2)

        plot!(Shape([Float64(x) for x in vcat(rounds, reverse(rounds))], 
                    [Float64(x) for x in vcat(lb, reverse(ub))]),
                    color=:lightgray,
                    line=nothing)

        hline!([cali[i]],
               color=:grey50,
               linestyle=:dash,
               linewidth=1.5)

        plot!(rounds, est,
              color=:black,
              linewidth=0.75)
    end

    for i = (length(cali)+1):subplot_tot
        plots[i] = plot(border = :none)
    end

    graph_dims = (gdim[2]*200, gdim[1]*115)
    plot_res = get_resulting_plot(plots, gdim, graph_dims, gord)

    savefig(plot_res, plotsdir(filename))
end


function plot_evolution_rmse(results, mod_set, folder)
    filename = "evolution_rmse_$(folder).pdf"

    cali = get_model_cali(mod_set)
    cons = get_model_cons(mod_set)
    gnam = get_model_gnam(mod_set)
    gdim = get_model_gdim(mod_set)
    gord = get_model_gord(mod_set)

    rsb, _ = results_best(results, cali)
    rounds = [sum(rsb[i][2]) for i in eachindex(rsb)]
    
    biases = []
    for i in eachindex(rsb)
        bias = []

        for j in eachindex(cali)
            bias = vcat(bias, abs.(rsb[i][1][j,:] .- cali[j]) ./ abs(cali[j]))
        end

        push!(biases, bias)
    end

    biases_mean = [mean(bias) for bias in biases]
    biases_sd = [1.960*std(bias)/sqrt(length(bias)) for bias in biases]
    #biases_sd = [std(bias) for bias in biases]

    biases_lb = biases_mean - biases_sd
    biases_ub = biases_mean + biases_sd

    plot = scatter(rounds, biases_mean,
                   legend=:none,
                   color=:black,
                   xlabel=L"number\ of\ moments",
                   ylabel=L"bias",
                   xticks=0:1:maximum(rounds),
                   ylim=(0.0, 1.0),
                   size=(1000,750),
                   left_margin=4mm,
                   bottom_margin=5mm)

    scatter!(rounds, biases_ub,
             marker=:hline,
             color=:black)
    scatter!(rounds, biases_lb,
             marker=:hline,
             color=:black)

    for i in eachindex(biases_sd)
        plot!([rounds[i],rounds[i]], [biases_lb[i], biases_ub[i]],
              color=:black,
              linestyle=:dash)
    end

    biases_argmin = argmin(biases_mean)

    efficient_ub = biases_mean .<= biases_ub[biases_argmin]
    efficient_rounds = rounds .>= length(cali)

    biases_efficient = argmax(efficient_ub .& efficient_rounds)

    vline!([rounds[biases_efficient]],
           color=:black,
           linestyle=:dash)

    hline!([biases_lb[biases_argmin], biases_ub[biases_argmin]],
           color=:black,
           linestyle=:dash)

    savefig(plot, plotsdir(filename))
end


##################
# PLOTS FOR TEXT #
##################

function plot_rmse_text(results, mod_set, folder; bench=[], bench_names=[], ylimits=nothing)
    filename = "text_rmse_$(folder).pdf"

    cali = get_model_cali(mod_set)

    results_rmse = get_rmse(results, cali)
    results_rounds = [sum(results[i][2]) for i in eachindex(results)]

    if isnothing(ylimits)
        plot = scatter(results_rounds, results_rmse,
                       markercolor=:white,
                       markeralpha=0.4,
                       legend=:none,
                       xlabel=L"\textit{number\ of\ moments}",
                       ylabel=L"\textit{RMSE}",
                       xguidefontsize=12,
                       yguidefontsize=12,
                       xticks=1:2:maximum(results_rounds),
                       xtickfontsize=10,
                       ytickfontsize=10,
                       xlims=(0.5,22.5),
                       size=(500,350),
                       left_margin=4mm,
                       bottom_margin=5mm)
    else
        plot = scatter(results_rounds, results_rmse,
                       markercolor=:white,
                       markeralpha=0.4,
                       legend=:none,
                       xlabel=L"\textit{number\ of\ moments}",
                       ylabel=L"\textit{RMSE}",
                       xguidefontsize=12,
                       yguidefontsize=12,
                       xticks=1:2:maximum(results_rounds),
                       xtickfontsize=10,
                       ytickfontsize=10,
                       xlims=(0.5,22.5),
                       ylims=ylimits,
                       size=(500,350),
                       left_margin=4mm,
                       bottom_margin=5mm)
    end

    rsb, best_idx = results_best(results, cali)
    rsb_annot = text.([" $(best_idx[i])" for i in 1:length(rsb)], 6, :left)
    rsb_rmse = get_rmse(rsb, cali)
    rsb_rounds = [sum(rsb[i][2]) for i in eachindex(rsb)]

    scatter!(rsb_rounds, rsb_rmse,
             series_annotations=rsb_annot,
             markercolor=:grey50)

    bench_annot = text.(["$i" for i in bench_names], 8, :bottom, :black)
    bench_rmse = get_rmse(bench, cali)
    bench_rounds = [sum(bench[i][2]) for i in eachindex(bench)]

    scatter!(bench_rounds, bench_rmse,
             series_annotations=bench_annot,
             markercolor=:black)

    savefig(plot, plotsdir(filename))
end


function plot_best_text(results, mod_set, folder; ylimits=nothing)
    filename = "text_best_$(folder).pdf"

    cali = get_model_cali(mod_set)

    results_rounds = [sum(results[i][2]) for i in eachindex(results)]

    rsb, best_idx = results_best(results, cali)
    rsb_annot = text.(["$(best_idx[i]) " for i in 1:length(rsb)], 6, :right)
    rsb_rmse = get_rmse(rsb, cali)
    rsb_rounds = [sum(rsb[i][2]) for i in eachindex(rsb)]

    if isnothing(ylimits)
        plot = scatter(rsb_rounds, rsb_rmse,
                       markercolor=:grey30,
                       legend=:none,
                       xlabel=L"\textit{number\ of\ moments}",
                       ylabel=L"\textit{RMSE_{\hat{\theta}}}",
                       xguidefontsize=12,
                       yguidefontsize=12,
                       xticks=1:2:maximum(rsb_rounds),
                       xtickfontsize=10,
                       ytickfontsize=10,
                       xlims=(0.5,22.5),
                       size=(500,200),
                       left_margin=4mm,
                       bottom_margin=5mm)
    else
        plot = scatter(rsb_rounds, rsb_rmse,
                       markercolor=:grey30,
                       legend=:none,
                       xlabel=L"\textit{number\ of\ moments}",
                       ylabel=L"\textit{RMSE_{\hat{\theta}}}",
                       xguidefontsize=12,
                       yguidefontsize=12,
                       xticks=1:2:maximum(rsb_rounds),
                       xtickfontsize=10,
                       ytickfontsize=10,
                       xlims=(0.5,22.5),
                       ylims=ylimits,
                       size=(500,200),
                       left_margin=4mm,
                       bottom_margin=5mm)
    end

    plot!(rsb_rounds, rsb_rmse,
          color=:grey30)

    savefig(plot, plotsdir(filename))
end


function plot_model_text(mod_set, obs, burn, type; seed=1, vline=nothing)
    filename = "$(mod_set["model"])_$(type).pdf"
    ylab = latexstring(type)

    cali = get_model_cali(mod_set)
    data = gen_series(mod_set, obs, burn, cali)

    Random.seed!(seed)

    p = plot(data,
             color=:black,
             legend=:none,
             xlab=L"time",
             ylab=ylab,
             size=(500,200))

    if !isnothing(vline)
        vline!([vline], color=:black, line=:dash)
    end

    savefig(p, plotsdir(filename))
end


function plot_evolution_bias_text(results, mod_set, folder, type, ylim=(0.0, 1.0))
    filename = "evolution_bias_$(folder).pdf"

    cali = get_model_cali(mod_set)
    cons = get_model_cons(mod_set)
    gnam = get_model_gnam(mod_set)
    gdim = get_model_gdim(mod_set)
    gord = get_model_gord(mod_set)

    rsb, _ = results_best(results, cali)

    if type == "bsw"
        rsb = reverse(rsb)
    end

    rounds = [sum(rsb[i][2]) for i in eachindex(rsb)]
    
    biases = []
    for i in eachindex(rsb)
        bias = []

        for j in eachindex(cali)
            bias = vcat(bias, abs.(rsb[i][1][j,:] .- cali[j]) ./ abs(cali[j]))
        end

        push!(biases, bias)
    end

    biases_mean = [mean(bias) for bias in biases]
    biases_sd = [1.960*std(bias)/sqrt(length(bias)) for bias in biases]
    #biases_sd = [std(bias) for bias in biases]

    biases_lb = biases_mean - biases_sd
    biases_ub = biases_mean + biases_sd

    biases_restrict = biases_ub - biases_lb .< ylim[2] - ylim[1]

    rounds = rounds[biases_restrict]
    rsb = rsb[biases_restrict]
    biases_mean = biases_mean[biases_restrict]
    biases_sd = biases_sd[biases_restrict]
    biases_lb = biases_lb[biases_restrict]
    biases_ub = biases_ub[biases_restrict]

    plot = scatter(rounds, biases_mean,
                   legend=:none,
                   markersize=3,
                   color=:black,
                   xlabel=L"\textit{number\ of\ moments}",
                   ylabel=L"\textit{B(\hat{\theta})}",
                   xguidefontsize=12,
                   yguidefontsize=12,
                   xticks=0:2:22,
                   xtickfontsize=10,
                   ytickfontsize=10,
                   xlims=(0.5,22.5),
                   ylim=ylim,
                   size=(500,350),
                   left_margin=4mm,
                   bottom_margin=5mm)

    scatter!(rounds, biases_ub,
             marker=:hline,
             color=:black)
    scatter!(rounds, biases_lb,
             marker=:hline,
             color=:black)

    for i in eachindex(biases_sd)
        plot!([rounds[i],rounds[i]], [biases_lb[i], biases_ub[i]],
              color=:black,
              linestyle=:dash)
    end

    biases_argmin = argmin(get_rmse(rsb, cali))

    hline!([biases_ub[biases_argmin]],
           color=:grey30,
           linestyle=:dash)

    efficient_ub = biases_mean .<= biases_ub[biases_argmin]
    efficient_rounds = rounds .>= length(cali)

    biases_efficient = argmax(efficient_ub .& efficient_rounds)

    vline!([rounds[biases_efficient]],
           color=:grey30,
           linestyle=:solid)

    savefig(plot, plotsdir(filename))
end


###############
# PLOTS UTILS #
###############

function get_resulting_plot(p, l, s, o)
    total = length(p)

    if total == 1
        result = plot(p[o[1]],
                      layout=l, size=s, right_margin=2mm)
    elseif total == 2
        result = plot(p[o[1]], p[o[2]],
                      layout=l, size=s)
    elseif total == 3
        result = plot(p[o[1]], p[o[2]], p[o[3]],
                      layout=l, size=s)
    elseif total == 6
        result = plot(p[o[1]], p[o[2]], p[o[3]], 
                      p[o[4]], p[o[5]], p[o[6]],
                      layout=l, size=s)
    elseif total == 9
        result = plot(p[o[1]], p[o[2]], p[o[3]], 
                      p[o[4]], p[o[5]], p[o[6]], 
                      p[o[7]], p[o[8]], p[o[9]],
                      layout=l, size=s)
    elseif total == 12
        result = plot(p[o[1]], p[o[2]], p[o[3]], 
                      p[o[4]], p[o[5]], p[o[6]], 
                      p[o[7]], p[o[8]], p[o[9]], 
                      p[o[10]], p[o[11]], p[o[12]],
                      layout=l, size=s)
    elseif total == 15
        result = plot(p[o[1]], p[o[2]], p[o[3]], 
                      p[o[4]], p[o[5]], p[o[6]], 
                      p[o[7]], p[o[8]], p[o[9]], 
                      p[o[10]], p[o[11]], p[o[12]], 
                      p[o[13]], p[o[14]], p[o[15]],
                      layout=l, size=s)
    elseif total == 18
        result = plot(p[o[1]], p[o[2]], p[o[3]], 
                      p[o[4]], p[o[5]], p[o[6]], 
                      p[o[7]], p[o[8]], p[o[9]], 
                      p[o[10]], p[o[11]], p[o[12]], 
                      p[o[13]], p[o[14]], p[o[15]], 
                      p[o[16]], p[o[17]], p[o[18]],
                      layout=l, size=s)
    else
        println("The generated plot size is undefined in get_resulting_plot().")
    end

    return result
end
