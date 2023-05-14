using Statistics
using JLD2
using PyPlot
pygui(true)

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


S_vals = [10, 20, 30, 40]

eq_out
eq_out[1]
eq_out[1][1].μ_M_ii
eq_out[1][1].ellipsis
eq_out[1][1].max_circle

plot(eq_out[1][1].μ_M_ii, eq_out[1][1].ellipsis, "o")
plot(eq_out[1][1].max_circle[1, :], eq_out[1][1].max_circle[2, :], "o")


