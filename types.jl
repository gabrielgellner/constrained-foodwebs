import Random: MersenneTwister
import LinearAlgebra: nullspace

# Utilities for State Variables and indexing
struct ForwardInverseVariable{T}
    values::Vector{T}
    idx_map::Dict{Tuple{Int, Int}, Int}
end

Base.getindex(fiv::ForwardInverseVariable, i::Int) = fiv.values[i]
Base.getindex(fiv::ForwardInverseVariable, i::Int, j::Int) = fiv.values[fiv.idx_map[(i, j)]]
function Base.setindex!(fiv::ForwardInverseVariable, value, i::Int)
    fiv.values[i] = value
    return
end
function Base.setindex!(fiv::ForwardInverseVariable, value, i::Int, j::Int)
    fiv.values[fiv.idx_map[(i, j)]] = value
    return
end

struct ForwardInverseVariableBound{T}
    bounds::Vector{T}
    idx_map::Dict{Tuple{Int, Int}, Int}
end
# Setup default constructor for bounds with single bound for all variables
function ForwardInverseVariableBound(idx_map, bound::T) where T
    return ForwardInverseVariableBound(fill(bound, length(idx_map)), idx_map)
end

Base.getindex(fivb::ForwardInverseVariableBound, i::Int) = fivb.bounds[i]
Base.getindex(fivb::ForwardInverseVariableBound, i::Int, j::Int) = fivb.bounds[fivb.idx_map[(i, j)]]
function Base.setindex!(fivb::ForwardInverseVariableBound, value, i::Int)
    fivb.bounds[i] = value
    return
end
function Base.setindex!(fivb::ForwardInverseVariableBound, value, i::Int, j::Int)
    fivb.bounds[fivb.idx_map[(i, j)]] = value
    return
end

# The Inverse to Forward mapping for interaction terms needed for the constructing an
# `ForwardInverseVariable`
function create_idx_map(net)
    idx_map = Dict{Tuple{Int, Int}, Int}()
    cur_i = 1
    for i in 1:size(net, 1)
        for j in 1:size(net, 2)
            if net[i, j] != 0
                idx_map[(i, j)] = cur_i
                cur_i += 1
            end
        end
    end
    return idx_map
end

function pmat_to_cnet(pmat::PredationMatrix, r)
    cnet = -pmat + pmat'
    for i in 1:pmat.S
        if r[i] > 0
            cnet[i, i] = -1
        elseif r[i] < 0
            cnet[i, i] = -1
        end
    end
    # # Assumption: all top species have no intraspecific terms -- i.e. have density independent mortality
    for s in top_species(pmat)
        cnet[s, s] = 0
    end
    return cnet
end

# Create inverse problem
function cnet_nz(cnet::Matrix)
    α_nz = [Tuple{Int, Int}[] for i in 1:size(cnet, 1)]
    for i in 1:size(cnet, 1)
        for j in 1:size(cnet, 2)
            if cnet[i, j] != 0
                push!(α_nz[i], (i, j))
            end
        end
    end
    return α_nz
end

function find_blocks(cnet, biomass)
    blocks = typeof(biomass)[]
    for i in 1:size(cnet, 1)
        push!(blocks, filter(x -> x != 0.0, cnet[i, :] .* biomass))
    end
    return blocks
end

function block_diagonal(α)
    block_size = length.(α)
    M = fill(0.0, length(α), sum(block_size))
    istart = 1
    for (i, block) in enumerate(α)
        iend = istart + block_size[i] - 1
        M[i, istart:iend] = block
        istart = iend + 1
    end
    return M
end

struct ConstrainedFoodWeb{T} <: AbstractMatrixConstraintProblem
    # Network
    S::Int
    cnet::Matrix{Int}
    # biomass and growth
    x::Vector{T}
    r::Vector{T}
    # equality coefficients
    E::Matrix{T}
    # inequality constraints
    # L * α >= f
    L::Matrix{T}
    f::Vector{T}
    # Inverse state variable (αi) to Forward matrix state variabel (Aij)
    idx_map::Dict{Tuple{Int, Int}, Int}
    # box constraints
    lower::ForwardInverseVariableBound{T}
    upper::ForwardInverseVariableBound{T}
end

function ConstrainedFoodWeb(pmat, x::Vector{T}, r::Vector{T}, lower::T, upper::T) where T
    @assert length(x) == length(r) == pmat.S
    
    cnet = pmat_to_cnet(pmat, r)
    idx_map = create_idx_map(cnet)
    
    E = block_diagonal(find_blocks(cnet, x))
    # Add in the default inequality constraint of a_prey <= a_pred
    f = fill(0.0, size(E, 2))
    L = fill(0.0, sum(pmat), size(E, 2))
    # The logic is we need a row of M for each predator / prey pair, and we want
    # a_ij - a_ji >= 0 (f)
    pred_count = 1
    for i in 1:pmat.S
        for j in 1:pmat.S
            if pmat[i, j] == 1
                L[pred_count, idx_map[(i, j)]] = 1
                L[pred_count, idx_map[(j, i)]] = -1
                pred_count += 1
            end
        end
    end
    
    # setup the bound types
    vlower = ForwardInverseVariableBound(idx_map, lower)
    vupper = ForwardInverseVariableBound(idx_map, upper)
    
    return ConstrainedFoodWeb(pmat.S, cnet, x, r, E, L, f, idx_map, vlower, vupper)
end

# Initialization
function ConstrainedFoodWeb(pmat; lb = 0.0, ub = 10.0, min_β = 10.0, δ = 0.0)
    x = set_β(pmat, δ, min_β)
    
    #TODO: likely this assumption about the `r` vector should be passed in, and then
    #      future calls to the model with set_r!(web) could change the `r` values inplace
    #      using the passed in method
    r = fill(-1.0, pmat.S)
    for xi in basal_species(pmat)
        r[xi] = 1.0
    end
    return ConstrainedFoodWeb(pmat, x, r, lb, ub)
end

# add a definition to generate a problem type for hit and run sampler
function HRSampler(conweb::ConstrainedFoodWeb; RNG = MersenneTwister())
    samp = HRSampler(conweb, nullspace(conweb.E), Vector{Vector{Float64}}(), fill(NaN, size(conweb.E, 2)), 0, 0, fill(NaN, size(conweb.E, 2)), 1e-15, 1e-15, false, RNG)
    warmup_sampler!(samp)
    return samp
end

function ACHRSampler(conweb::ConstrainedFoodWeb; RNG = MersenneTwister())
    samp = ACHRSampler(conweb, nullspace(conweb.E), Vector{Vector{Float64}}(), fill(NaN, size(conweb.E, 2)), 0, 0, fill(NaN, size(conweb.E, 2)), 1e-15, 1e-15, false, RNG)
    warmup_sampler!(samp)
    return samp
end

function α_mat(conweb::ConstrainedFoodWeb, αs)
    A = fill(0.0, conweb.S, conweb.S)
    for ((i, j), α_i) in conweb.idx_map
        A[i, j] = conweb.cnet[i, j] * αs[α_i]
    end
    return A
end

cmat(conweb::ConstrainedFoodWeb, samp) = α_mat(conweb, samp) .* conweb.x