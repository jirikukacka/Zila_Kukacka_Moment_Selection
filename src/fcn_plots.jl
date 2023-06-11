using JLD, Plots, Plots.PlotMeasures, StatsPlots

include("fcn_watson.jl")
include("fcn_metric.jl")
include("fcn_utils.jl")
include("fcn_data.jl")
include("list_models.jl")


function plot_rmse(results, mod_set, folder; bench=[], bench_names=[], ylimits=nothing, show_plot=false)
    filename = "rmse_$(folder).pdf"

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

    if show_plot
        display(plot)
    end

    savefig(plot, plotsdir(filename))
end

function plot_evolution_bias(results, mod_set, folder, type, ylim=(0.0, 1.0); show_plot=false)
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

    if show_plot
        display(plot)
    end

    savefig(plot, plotsdir(filename))
end

function plot_evolution_par(results, mod_set, folder; bench=[], bench_names=[], ylims=nothing, show_plot=false)
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
    
    if show_plot
        display(plot_res)
    end
    
    savefig(plot_res, plotsdir(filename))
end

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
