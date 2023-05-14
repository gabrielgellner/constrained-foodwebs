using LinearAlgebra
using Distributions
using Statistics
using PyPlot
pygui(true)

include("foodwebs.jl")
include("biomass_tools.jl")

function rand_stab(pmat, dist; f = 0.8)
    M = fill(0.0, pmat.S, pmat.S)
    top_idx = top_species(pmat)
    for i in 1:pmat.S
        for j in 1:pmat.S
            if i == j
                if i in top_idx
                    M[i, j] == 0.0
                else
                    M[i, j] = rand(Uniform(-1.1, -0.8))
                end
            elseif pmat[i, j] == 1
                M[i, j] = -rand(dist)
                M[j, i] = f * rand(dist)
            end
        end
    end

    return λ_stability(M)
end

function div_stab(S_vals, C, dist)
    S_stab = fill(0.0, 3, length(S_vals))
    n_net = 12 * 42
    for (i_S, S) in enumerate(S_vals)
        @show S
        net_stab = fill(0.0, n_net)
        for i_net in 1:n_net
            pmat = generalized_cascade_network(S, C)

            n_samp = 100
            stab = fill(0.0, n_samp)
            for i in 1:n_samp
                stab[i] = rand_stab(pmat, dist)
            end

            net_stab[i_net] = mean(stab)
        end
        S_stab[1, i_S] = mean(net_stab)
        S_stab[2, i_S] = mean(net_stab) - std(net_stab)
        S_stab[3, i_S] = mean(net_stab) + std(net_stab)
    end
    return S_stab
end

# Full space is over (0, 1)
dist = Uniform(0, 5.0)
C = 0.25
S_vals = 10:2:40

S_stab = div_stab(S_vals, C, dist)

# Equal Scale
let
    figure()
    plot(S_vals, S_stab[1, :])
    plot(S_vals, S_stab[2, :], "k")
    plot(S_vals, S_stab[3, :], "k")

    axhline(y = 0, linestyle = "--", color = "gray")
    xlabel("Diversity: S")
    ylabel("Stability: λ1")
    ylim(-1, 5)
end

# Zoom in
let
    figure()
    plot(S_vals, S_stab[1, :])
    plot(S_vals, S_stab[2, :], "k")
    plot(S_vals, S_stab[3, :], "k")

    axhline(y = 0, linestyle = "--", color = "gray")
    xlabel("Diversity: S")
    ylabel("Stability: λ1")

end
