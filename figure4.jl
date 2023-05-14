using Statistics
using JLD2
using PyPlot
pygui(true)

function make_feed_dist(data; M_ij = true)
    # First find the maximum # of feeding links in a set
    n_feed_max = maximum(maximum(out.link_dist[:, 1, 1]) for out in data)
    feed_dist = Dict(1:n_feed_max .=> [Float64[] for i in 1:round(Int, n_feed_max)])

    for i in eachindex(data)
        for (n, n_feed) in enumerate(data[i].link_dist[:, 1, 1])
            if n_feed === NaN || n_feed == 0.0
                continue
            end
            if M_ij
                push!(feed_dist[n_feed], data[i].link_dist[n, 2, 1])
            else
                push!(feed_dist[n_feed], data[i].link_dist[n, 3, 1])
            end
        end
    end

    return feed_dist
end

# Data locations
uncon_dir = "data/uncon"
uncon_files =  readdir(uncon_dir)

eqcon_dir = "data/eqcon"
eqcon_files = readdir(eqcon_dir)

ieqcon_dir = "data/ieqcon"
ieqcon_files = readdir(ieqcon_dir)

## Unconstrained
data = load(joinpath(uncon_dir, uncon_files[1])) 
S_vals = data["S_vals"]
uncon_out = data["out"]

## Constrained
data = load(joinpath(eqcon_dir, eqcon_files[1]))
eq_out = data["out"]

## Constrained + Energetic
data = load(joinpath(ieqcon_dir, ieqcon_files[1]))
ieq_out = data["out"]


function plot_feed_dist(feed_dist; minmax=true, var_bounds=false, show_x_label=true)
    summary = fill(0.0, length(feed_dist), 5)
    for (i, (k, v)) in enumerate(feed_dist)
        if length(v) == 0
            continue
        end
        summary[i, 1] = k
        summary[i, 2] = minimum(v)
        summary[i, 3] = maximum(v)
        summary[i, 4] = mean(v)
        summary[i, 5] = std(v)
    end
    
    if minmax
        plot(summary[:, 1], summary[:, 2], "ko")
        plot(summary[:, 1], summary[:, 3], "ko")
    end
    if var_bounds
        plot(summary[:, 1], summary[:, 4] .- summary[:, 5], "ko")
        plot(summary[:, 1], summary[:, 4] .+ summary[:, 5], "ko")
    end
    plot(summary[:, 1], summary[:, 4], "o")
    if show_x_label
        xlabel("Number of Prey Links")
    end
    ylabel("Average Mean IS")
    
    return
end



let
    matplotlib.rc("font", size = 12)
    figure(figsize=(8, 10))
    # Uncon plot
    subplot(3, 2, 1)
    plot_feed_dist(make_feed_dist(uncon_out[4]); minmax=false, show_x_label=false)
    title("Negative Impact of Pred on Prey")
    xlim(0.5, 16.5)
    ylim(1, 3)
    
    subplot(3, 2, 2)
    plot_feed_dist(make_feed_dist(uncon_out[4], M_ij=false); minmax=false, show_x_label=false)
    title("Positive Impact of Prey on Pred")
    xlim(0.5, 16.5)
    ylim(1, 3)

    # Feasible plot
    subplot(3, 2, 3)
    plot_feed_dist(make_feed_dist(eq_out[4]); minmax=false, show_x_label=false)
    xlim(0.5, 16.5)
    ylim(0, 47)
    
    subplot(3, 2, 4)
    plot_feed_dist(make_feed_dist(eq_out[4], M_ij=false); minmax=false, show_x_label=false)
    xlim(0.5, 16.5)
    ylim(0.46)

    # Feasible + Energ
    subplot(3, 2, 5)
    plot_feed_dist(make_feed_dist(ieq_out[4]); minmax=false)
    xlim(0.5, 16.5)
    ylim(0, 1.4)
    
    subplot(3, 2, 6)
    plot_feed_dist(make_feed_dist(ieq_out[4], M_ij=false); minmax=false)
    xlim(0.5, 16.5)
    ylim(0, 1.4)

    tight_layout()

    savefig("figs/figure4_only_mean.png", dpi=300)
end

