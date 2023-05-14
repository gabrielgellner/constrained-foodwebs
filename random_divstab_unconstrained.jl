using LinearAlgebra
using Distributions
using Statistics
using PyPlot
pygui(true)

include("foodwebs.jl")
include("biomass_tools.jl")

function rand_stab(pmat, β, dist)
    A = fill(0.0, pmat.S, pmat.S)
    top_idx = top_species(pmat)
    for i in 1:pmat.S
        for j in 1:pmat.S
            if i == j
                if i in top_idx
                    A[i, j] == 0.0
                else
                    A[i, j] = -rand(dist)
                end
            elseif pmat[i, j] == 1
                A[i, j] = -rand(dist)
                A[j, i] = rand(dist)
            end
        end
    end

    M = A .* β
    return λ_stability(M)
end

function div_stab(S_vals, C, dist)
    S_stab = fill(0.0, 3, length(S_vals))
    for (i_S, S) in enumerate(S_vals)
        n_net = 12 * 42
        net_stab = fill(0.0, n_net)
        for i_net in 1:n_net
            pmat = generalized_cascade_network(S, C)
            min_β = 10.0

            n_samp = 100
            stab = fill(0.0, n_samp)
            for i in 1:n_samp
                δ = rand(Uniform(-5.0, 5.0))
                β = set_β(pmat, δ, min_β)

                stab[i] = rand_stab(pmat, β, dist)
            end

            net_stab[i_net] = mean(stab)
        end
        S_stab[1, i_S] = mean(net_stab)
        S_stab[2, i_S] = mean(net_stab) - std(net_stab)
        S_stab[3, i_S] = mean(net_stab) + std(net_stab)
    end
    return S_stab
end

dist = Uniform(0, 1.0)
C = 0.25
S_vals = 10:40

S_stab = div_stab(S_vals, C, dist)

let
    plot(S_vals, S_stab[1, :])
    plot(S_vals, S_stab[2, :], "k")
    plot(S_vals, S_stab[3, :], "k")

    axhline(y = 0, linestyle = "--", color = "gray")
    xlabel("Diversity: S")
    ylabel("Stability: λ1")
end
