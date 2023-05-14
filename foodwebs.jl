using LinearAlgebra
using Distributions
using RecursiveArrayTools
import Graphs: DiGraph, is_connected

λ_stability(M) = maximum(real.(eigvals(M)))

struct PredationMatrix <: AbstractArray{Int, 2}
    links::Matrix{Int}
    S::Int
    C::Float64
end

function PredationMatrix(links::Matrix{Int})
    S1, S2 = size(links)
    if S1 != S2
        error("PredationMatrix must be square, given dimensions ($S1, S2)")
    end
    # if there is canibals the connectance measure might be off
    return PredationMatrix(links, S1, 2 * sum(links) / S1 ^ 2)
end

Base.size(pmat::PredationMatrix) = size(pmat.links)
Base.IndexStyle(::Type{<:PredationMatrix}) = IndexCartesian()
Base.getindex(pmat::PredationMatrix, I...) = getindex(pmat.links, I...)
Base.strides(pmat::PredationMatrix) = strides(pmat.links)

"""
    generalized_cascade_network(S, C)::PredationMatrix

Generate a `S`×`S` adjancency matrix, with connectance `C`, following the network structure
described by:
Stouffer, D. B., Camacho, J., Guimerà, R. & Ng, C. Quantitative patterns in the structure of
model and empirical food webs. Ecology 86, 1301–1311 (2005).
"""
function generalized_cascade_network(S::Int, C::Float64; maxiter = 1000)
    for iter in 1:maxiter
        β = 1 / C - 1
        adj = fill(0, S, S)
        for i = 1:S
            for j = (i + 1):S
                if rand() < rand(Beta(1, β))
                    adj[i, j] = 1
                end
            end
        end
        if is_connected(DiGraph(adj))
            return PredationMatrix(adj)
        end
    end
    error("Could not find a connected network for ($S, $C) in $maxiter iteractions")
end

connectance(adj) = sum(adj) / size(adj, 2) ^ 2
connectance(adj::PredationMatrix) = 2 * sum(adj) / adj.S ^ 2

function basal_species(pmat::PredationMatrix)
    col_sums = sum(pmat, dims = 1)
    return [i for i in 1:length(col_sums) if col_sums[i] == 0]
end

function top_species(pmat::PredationMatrix)
    row_sums = sum(pmat, dims = 2)
    return [i for i in 1:length(row_sums) if row_sums[i] == 0]
end

function trophic_levels(adj)
    # the matrix comes in "predator" form but the matrix algebra requires it to be in
    # "prey" form, also needs to be float since we use transition probabilites
    pmat = collect(1.0 * adj')
    for i in 1:size(pmat, 1)
        pmat[i, i] = 0
    end

    for i in 1:size(pmat, 1)
        total_prey = sum(pmat[i, :])
        if total_prey != 0
            pmat[i, :] ./= total_prey
        end
    end

    return pinv(I - pmat) * fill(1.0, size(pmat, 1))
end

"""
    offdiag_pairs(mat::Matrix)

return an array of the pairs of the offdiagonal pairs of the square matrix `mat`

"""
function offdiag_pairs(mat::Matrix)
    S = size(mat, 1)
    @assert S == size(mat, 2) # must be square
    pairs = fill(zero(eltype(mat)), Int(S * (S - 1) / 2), 2)
    current_row = 1
    for i = 1:S
        for j = i:S
            if i != j
                pairs[current_row, :] = [mat[i, j], mat[j, i]]
                current_row += 1
            end
        end
    end
    return pairs
end

function offdiag_pairs_nz(M::Matrix)
    out = []
    for i in 1:(size(M, 1) - 1)
        for j in (i + 1):size(M, 2)
            if M[i, j] != 0
                push!(out, [M[i, j], M[j, i]])
            end
        end
    end
    return VectorOfArray(out)
end