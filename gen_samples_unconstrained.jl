using Random
using JLD2

include("biomass_tools.jl")
include("foodwebs.jl")

function rand_M(pmat, β, dist; sym = true, f = 0.8)
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
                if sym
                    A[j, i] = rand(dist)
                else
                    A[j, i] = f * rand(dist)
                end
            end
        end
    end

    # M matrix
    return A .* β
end

function calc_link_dist(S, M)
    n_feed_link = fill(0, S)
    μ_IS = fill(0.0, S)

    # take the mean of the columns with negative elements, but not diagonals
    for j in 1:S
        for i in 1:S
            if i == j
                continue
            end
            if M[i, j] < 0
                n_feed_link[j] += 1
                μ_IS[j] +=  abs(M[i, j])
            end
        end
    end
    μ_IS = μ_IS ./ n_feed_link
    return [n_feed_link μ_IS]
end

function rand_summary(S, n_samp; sym = true, ub = 1.0)
    # replicate the output of the full sampler for ease of use
    out = fill(0.0, 12, n_samp)
    out = (
        μ_M_ij = fill(NaN, n_samp),
        σ_M_ij = fill(NaN, n_samp),
        λ1 = fill(NaN, n_samp),
        μ_M_ii = fill(NaN, n_samp),
        σ_M_ii = fill(NaN, n_samp),
        # Note every chain has a SINGLE network associated with it
        #conweb = deepcopy(c.conweb)
        # So we have 
        # will have S entries with a (#link, mean(IS-))
        link_dist = fill(NaN, S, 2, n_samp)
    )

    pmat = generalized_cascade_network(S, 0.25)
    δ = rand(Uniform(-5, 5))
    min_β = 10.0
    β = set_β(pmat, δ, min_β)
    for i in 1:n_samp

        M = rand_M(pmat, β, Uniform(0, ub); sym = sym)
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
        out.link_dist[:, :, i] = calc_link_dist(S, M)
    end

    return out
end

uncon_sym = [[rand_summary(S, 1000) for i in 1:500] for S in [20, 30, 40, 50, 60]]
uncon = [[rand_summary(S, 500; sym = false) for i in 1:100] for S in [10, 20, 30, 40]]

function main(S_vals, mul_samp; sym = true, ub = 1.0)
    out = []
    for (i, S) in enumerate(S_vals)
        @show S
        # If we make the n_samp a multiple of Threads.nthreads() (number of processsors on computer, likely be the most efficient)
        push!(out, rand_summary(S, Threads.nthreads() * mul_samp; sym = sym))
    end
    return out
end
    
let
    S_vals = [10, 20, 30, 40]
    n_replicate = 1
    out = main(S_vals, n_replicate; sym = false, ub = 5.0)

    save("data/uncon/con_slice" * randstring() * ".jld2", "out", out, "S_vals", S_vals)
end