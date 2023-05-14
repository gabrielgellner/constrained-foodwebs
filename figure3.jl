using Statistics
using JLD2
using PyPlot
using Seaborn
pygui(true)

function slice_plot(out, slice_i, label; loc = "upper right")
    figure()
    for (i, S_out) in enumerate(out)
        S = S_vals[i]
        kdeplot([mean(getfield(o, slice_i)) for o in S_out], label = "S = $S", fill = true) #, kde_kws = Dict("shade" => true, "linewidth" => 2))
    end
    xlabel(label)
    legend(loc = loc)
end

# Data locations
uncon_dir = "data/uncon"
uncon_files =  readdir(uncon_dir)

eqcon_dir = "data/eqcon"
eqcon_files = readdir(eqcon_dir)

ieqcon_dir = "data/ieqcon"
ieqcon_files = readdir(ieqcon_dir)


# Figure 3a: IS Structure
## Unconstrained
### Density plot of Mij averages
data = load(joinpath(uncon_dir, uncon_files[1]))
S_vals = data["S_vals"]
uncon_out = data["out"]
slice_plot(uncon_out, :μ_M_ij, "Mean Mij")
title("Unconstrained")
savefig("figs/figure3a1.png")

## Feasible Constraints
### Density plot of Mij averages
data = load(joinpath(eqcon_dir, eqcon_files[1]))
S_vals = data["S_vals"]
eq_out = data["out"]
slice_plot(eq_out, :μ_M_ij, "Mean Mij")
title("Feasible Constraints")
savefig("figs/figure3a2.png")

## Feasible + Energetic Constraints
### Density plot of Mij averages
data = load(joinpath(ieqcon_dir, ieqcon_files[1]))
S_vals = data["S_vals"]
ieq_out = data["out"]
slice_plot(ieq_out, :μ_M_ij, "Mean Mij")
xlim(-0.5, 5)
title("Feasible + Energetic Constraints")
savefig("figs/figure3a3.png")


# Figure 3b: Diagonal Structure Mean
## Unconstrained
slice_plot(uncon_out, :μ_M_ii, "Mean Mii")
title("Unconstrained")
savefig("figs/figure3b1.png", dpi = 300)

## Feasible Constraints
slice_plot(eq_out, :μ_M_ii, "Mean Mii"; loc = "upper left")
title("Feasible Constraints")
savefig("figs/figure3b2.png")

## Energetic Constraints
slice_plot(ieq_out, :μ_M_ii, "Mean Mii"; loc = "upper left")
xlim(-1, 0.1)
title("Feasible + Energetic Constraints")
savefig("figs/figure3b3.png")


# Figure 3c: Diagonal Structure Variance
## Unconstrained
slice_plot(uncon_out, :σ_M_ii, "Variance Mii"; loc = "upper left")
xlim(0.084, 0.087)
title("Unconstrained")
savefig("figs/figure3c1.png")

## Feasible Constraints
slice_plot(eq_out, :σ_M_ii, "Variance Mii")
title("Feasible Constraints")
savefig("figs/figure3c2.png")

## Energetic Constraints
slice_plot(ieq_out, :σ_M_ii, "Variance Mii")
xlim(-0.1, 2)
title("Feasible + Energetic Constraints")
savefig("figs/figure3c3.png")
