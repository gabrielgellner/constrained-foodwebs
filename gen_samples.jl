using Random
using JLD2

include("sampler.jl")
include("biomass_tools.jl")
include("foodwebs.jl")
include("types.jl")

const mts = Random.MersenneTwister.(i for i in 1:Threads.nthreads())

function diversity_samples(pmat; max_attempt = 1000,  energetic = true, ub = 5.0, δ_min = -5.0, δ_max = 5.0, min_β = 10.0, basal_max = 5.0, nonbasal_max = 1.0)
    mt = mts[Threads.threadid()]
    δ = rand(mt, Uniform(δ_min, δ_max))

    conweb = ConstrainedFoodWeb(pmat; δ = δ, min_β = min_β, lb = 0.0, ub = ub)
    R_index = basal_species(pmat)
    if !energetic
        conweb.L[:] .= 0.0
    end

    sampler = ACHRSampler(conweb; RNG = mt)
    for attempt in 1:max_attempt
        for i in 1:pmat.S
            if i in R_index
                conweb.r[i] = rand(mt, Uniform(0, basal_max))
            else
                conweb.r[i] = -rand(mt, Uniform(0, nonbasal_max))
            end
        end
        ## Setup the sampler
        sampler = ACHRSampler(conweb; RNG = mt)
        if sampler.feasible
            #@show attempt
            break
        end
    end
    chain = shuffle!(hr_sample!(sampler, 1_000_000, thinning = 10))
    return (feasible = true, conweb = conweb, αs = chain)
end

function calc_link_dist(S, M)
    n_feed_link = fill(0, S)
    μ_IS = fill(0.0, S)
    μ_ISp = fill(0.0, S)

    # take the mean of the columns with negative elements, but not diagonals
    for j in 1:S
        for i in 1:S
            if i == j
                continue
            end
            if M[i, j] < 0
                n_feed_link[j] += 1
                μ_IS[j] +=  abs(M[i, j])
                μ_ISp[j] +=  abs(M[j, i])
            end
        end
    end
    μ_IS = μ_IS ./ n_feed_link
    μ_ISp = μ_ISp ./ n_feed_link
    return [n_feed_link μ_IS μ_ISp]
end

function row_sums(A)
    n_row, n_col = size(A)
    centers = zeros(eltype(A), n_row)
    rsums = zeros(eltype(A), n_row)
    for i in 1:n_row
        for j in 1:n_col
            if i == j
                centers[i] = A[i, j]
                continue
            end
            rsums[i] += abs(A[i, j])
        end
    end
    return (centers, rsums)
end

function max_circle(A)
    centers, rsums = row_sums(A)
    i_max = argmax(centers .+ rsums)
    μ_circle = mean(rsums)
    return (centers[i_max], rsums[i_max], μ_circle)
end

function ellipsis_measure(S, pairs)
    σ2x, ρ21, σ2y = cov(pairs, dims = 1, corrected = false)[[1, 2, 4]]
    σ2 = sqrt(σ2x * σ2y)
    ρ = ρ21 / σ2
    return sqrt(S * σ2) * (1 + ρ)
end

function summary_chain(c)
    #TODO: turn this into a struct
    out = (
        μ_M_ij = fill(NaN, length(c.αs)),
        σ_M_ij = fill(NaN, length(c.αs)),
        λ1 = fill(NaN, length(c.αs)),
        μ_M_ii = fill(NaN, length(c.αs)),
        σ_M_ii = fill(NaN, length(c.αs)),
        link_dist = fill(NaN, c.conweb.S, 3, length(c.αs)),
        # Gerschgorin measures
        max_circle_center = fill(NaN, length(c.αs)),
        max_circle_width = fill(NaN, length(c.αs)),
        μ_circle_width = fill(NaN, length(c.αs)),
        # Elipsis
        ellipsis = fill(NaN, length(c.αs)) 
    )

    for (i, α) in enumerate(c.αs)
        M = cmat(c.conweb, α)
        pairs = offdiag_pairs_nz(M)

        # μ_M_ij[i]
        out.μ_M_ij[i] = mean(abs.(pairs[1, :]))
        # μ_M_ji[i]
        out.σ_M_ij[i] = std(abs.(pairs[1, :]))
        #out[7, i] = λ_stability(M)
        out.λ1[i] = λ_stability(M)
        # Mean and Var of Diagonal entries
        #out[11, i] = mean(diag(M))
        out.μ_M_ii[i] = mean(diag(M))
        #out[12, i] = var(diag(M))
        out.σ_M_ii[i] = std(diag(M))

        # link dist
        out.link_dist[:, :, i] = calc_link_dist(c.conweb.S, M)

        # gershgorin time
        center, width, μ_circle_width = max_circle(M)
        out.max_circle_center[i] = center
        out.max_circle_width[i] = width
        out.μ_circle_width[i] = μ_circle_width

        # ellipsis
        out.ellipsis[i] = ellipsis_measure(c.conweb.S, pairs)
    end

    return out
end

function summary_stab(S, C, n_samp; energetic = true, ub = 5.0)
    out = []
    Threads.@threads for i in 1:n_samp
        chain = diversity_samples(generalized_cascade_network(S, C); energetic = energetic, ub = ub)
        
        if sum(isnan.(chain.αs)) == 0
           push!(out, summary_chain(chain))
        end
    end
    return out
end


# Diversity sample
function main(S_vals, C, mul_samp; energetic = true, ub = 1.0)
    out = []
    for (i, S) in enumerate(S_vals)
        @show S
        # If we make the n_samp a multiple of Threads.nthreads() (number of processsors on computer, likely be the most efficient)
        push!(out, summary_stab(S, C, Threads.nthreads() * mul_samp; energetic, ub))
    end
    return out
end
    
let
    C = 0.25
    S_vals = [10, 20, 30, 40]
    # Really this is the number of different (network, r, biomass) triples we are generating the samples for at each S
    # though multiplied by # of cores ... so in this case 12*, so for 100 we have 1200 different (network, biomass, r-vector) triples we explore
    n_replicate = 100
    out = main(S_vals, C, n_replicate; energetic = false, ub = 5.0)

    save("data/eqcon/con_slice" * randstring() * ".jld2", "out", out, "S_vals", S_vals)
end

let
    C = 0.25
    S_vals = [10, 20, 30, 40]
    n_replicate = 100
    out = main(S_vals, C, n_replicate; energetic = true, ub = 5.0)

    save("data/ieqcon/con_slice_energetic" * randstring() * ".jld2", "out", out, "S_vals", S_vals)
end  