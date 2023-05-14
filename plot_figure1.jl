using CSV
using DataFrames
using PyPlot
pygui(true)

data = CSV.read("figure_01.csv", DataFrame)

S_vals = 10:40

let
    plot(S_vals, data[:, 1])
    plot(S_vals, data[:, 2], "k")
    plot(S_vals, data[:, 3], "k")
    axhline(y = 0, linestyle = "--", color = "gray")
    xlabel("Diversity: S")
    ylabel("Stability: λ1")
    ylim(-1, 5)
end