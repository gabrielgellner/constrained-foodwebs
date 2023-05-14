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
eqcon_dir = "data/eqcon"
eqcon_files = readdir(eqcon_dir)

ieqcon_dir = "data/ieqcon"
ieqcon_files = readdir(ieqcon_dir)

## Feasible Constraints
### Density plot of Mij averages
data = load(joinpath(eqcon_dir, eqcon_files[1]))
S_vals = data["S_vals"]
eq_out = data["out"]

slice_plot(eq_out, :μ_M_ii, "Mean Mii")
title("Feasible Constraints")
savefig("figs/figure5a1.png") 

slice_plot(eq_out, :ellipsis, "Ellipsis")
title("Feasible Constraints")
savefig("figs/figure5a2.png") 

slice_plot(eq_out, :max_circle_center, "Max Circle Center"; loc = "upper")
#xlim(-0.5, 1)
title("Feasible Constraints")
savefig("figs/figure5a3.png") 

slice_plot(eq_out, :max_circle_width, "Max Circle Width")
title("Feasible Constraints")
savefig("figs/figure5a4.png") 

slice_plot(eq_out, :μ_circle_width, "Mean Circle Width")
title("Feasible Constraints")
savefig("figs/figure5a5.png") 

let
    out = eq_out
    figure()
    for (i, S_out) in enumerate(out)
        S = S_vals[i]
        kdeplot([mean(getfield(o, :μ_circle_width) .- getfield(o, :μ_M_ii)) for o in S_out], label = "S = $S", fill = true)
    end
    title("Feasible Constraints")
    xlabel("Mean Gershgorin Circle")
    legend()
    savefig("figs/figure5a6.png")     
end

## Feasible + Energetic Constraints
### Density plot of Mij averages
data = load(joinpath(ieqcon_dir, ieqcon_files[1]))
S_vals = data["S_vals"]
ieq_out = data["out"]

slice_plot(ieq_out, :μ_M_ii, "Mean Mii")
xlim(-1.0, 1)
title("Feasible + Energetic Constraints")
savefig("figs/figure5b1.png") 


slice_plot(ieq_out, :ellipsis, "Ellipsis")
title("Feasible + Energetic Constraints")
savefig("figs/figure5b2.png") 


slice_plot(ieq_out, :max_circle_center, "Max Circle Center")
xlim(-0.5, 1)
title("Feasible + Energetic Constraints")
savefig("figs/figure5b3.png") 


slice_plot(ieq_out, :max_circle_width, "Max Circle Width")
title("Feasible + Energetic Constraints")
savefig("figs/figure5b4.png") 

slice_plot(ieq_out, :μ_circle_width, "Mean Circle Width")
title("Feasible + Energetic Constraints")
savefig("figs/figure5b5.png") 

let
    out = ieq_out
    figure()
    for (i, S_out) in enumerate(out)
        S = S_vals[i]
        kdeplot([mean(getfield(o, :μ_circle_width) .- getfield(o, :μ_M_ii)) for o in S_out], label = "S = $S", fill = true)
    end
    xlabel("Mean Gershgorin Circle")
    legend()
    title("Feasible + Energetic Constraints")
    savefig("figs/figure5b6.png")     
end