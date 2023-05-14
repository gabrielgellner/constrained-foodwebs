using Random
using PyPlot
pygui(true)

include("foodwebs.jl")
include("sampler.jl")
include("biomass_tools.jl")
include("types.jl")

const mts = Random.MersenneTwister.(i for i in 1:Threads.nthreads())

function diversity_stab(pmat; max_attempt = 1000,  energetic = true, ub = 5.0)
    mt = mts[Threads.threadid()]
    δ = rand(mt, Uniform(-5.0, 5.0))
    min_β = 10.0

    conweb = ConstrainedFoodWeb(pmat; δ = δ, min_β = min_β, lb = 0.0, ub = ub)
    R_index = basal_species(pmat)
    if !energetic
        conweb.L[:] .= 0.0
    end
    #NOTE: this is done to avoid error of not being defined inside for loop -- think of better way
    sampler = ACHRSampler(conweb; RNG = mt)
    for attempt in 1:max_attempt
        for i in 1:pmat.S
            if i in R_index
                conweb.r[i] = rand(mt, Uniform(0, 5))
            else
                conweb.r[i] = -rand(mt, Uniform(0, 1))
            end
        end
        ## Setup the sampler
        sampler = ACHRSampler(conweb; RNG = mt)
        if sampler.feasible
            #@show attempt
            break
        end
    end
    chain = shuffle!(hr_sample!(sampler, 100_000, thinning = 10))
    
    return (feasible = true, conweb = conweb, αs = chain)
end

function summary_chain(c)
    out = fill(0.0, 3, length(c.αs))

    for (i, α) in enumerate(c.αs)
        M = cmat(c.conweb, α)
        out[1, i] = λ_stability(M)
    end

    return out
end

function summary_stab(S; n_samp = 10, energetic = true, ub = 5.0)
    out = []
    Threads.@threads for i in 1:n_samp
        chain = diversity_stab(generalized_cascade_network(S, 0.25); energetic = energetic, ub = ub)
        if sum(isnan.(chain.αs)) == 0
            push!(out, summary_chain(chain))
        end
    end
    return out
end

# Inner and outer variance
function main(S_vals)
    out = fill(NaN, 5, length(S_vals))
    for (i, S) in enumerate(S_vals)
        @show S
        chains = summary_stab(S; n_samp = Threads.nthreads() * 42, energetic = true, ub = 1.0)
        μ = mean(mean(c[1, :]) for c in chains)
        # So this gives the variance of the mean stability of each "network starting point" versus above
        σ = std(mean(c[1, :]) for c in chains)
        # which is just the variance for a given set of samples
        σ_inner = mean(std(c[1, :]) for c in chains)
        out[1, i] = μ
        out[2, i] = μ - σ
        out[3, i] = μ + σ
        out[4, i] = μ - σ_inner
        out[5, i] = μ + σ_inner
    end
    return out
end

S_vals = 10:1:20
out = main(S_vals)

let
    plot(S_vals, out[1, :])
    plot(S_vals, out[2, :], "k")
    plot(S_vals, out[3, :], "k")
    axhline(y = 0, linestyle = "--", color = "gray")
    xlabel("Diversity: S")
    ylabel("Stability: λ1")
    ylim(-1, 5)
end

using CSV
using Tables
using DataFrames

CSV.write("figure_01_eng.csv", Tables.table(out'))

data = CSV.read("figure_01.csv")
size(data)
plot(S_vals, data[1, :])