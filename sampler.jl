using LinearAlgebra
using Random
using Distributions
using JuMP
using Gurobi
#using Clp
using RecursiveArrayTools

abstract type AbstractMatrixConstraintProblem end
# This needs fields for E * α + r = 0:
# biomass (x) and growth (r)
# - x::Vector{Float64}
# - r::Vector{Float64}
# box constraints
# - lower::Vector{Float64}
# - upper::Vector{Float64}
# equality coefficients
# - E::Matrix{Float64}
# inequality M α >= f
# - L::Matrix{Float64}
# - f::Vector{Float64}

mutable struct HRSampler{T <: AbstractMatrixConstraintProblem}
    # Problem information
    prob::T
    Z::Matrix{Float64} # nullspace of prob.E
    # Chain initialization
    warmup::Vector{Vector{Float64}}
    center::Vector{Float64}
    n_warmup::Int
    n_samples::Int # how many samples have already been generated, used in center calculations
    prev::Vector{Float64}
    # Numerical checks
    feasibility_tol::Float64 # 1e-15 is a good start
    bounds_tol::Float64
    feasible::Bool
    RNG::Random.MersenneTwister
end

mutable struct ACHRSampler{T <: AbstractMatrixConstraintProblem}
    # Problem information
    prob::T
    Z::Matrix{Float64} # nullspace of prob.E
    # Chain initialization
    warmup::Vector{Vector{Float64}}
    center::Vector{Float64}
    n_warmup::Int
    n_samples::Int # how many samples have already been generated, used in center calculations
    prev::Vector{Float64}
    # Numerical checks
    feasibility_tol::Float64 # 1e-15 is a good start
    bounds_tol::Float64
    feasible::Bool
    RNG::Random.MersenneTwister
end

is_feasible(prob, α; atol = 1e-15) = all(isapprox.(prob.E * α + prob.r, 0.0, atol = atol))

# build a constant environment to not restart the gurobi environment for each call
const GRB_ENV = Gurobi.Env(output_flag=0)

function warmup_sampler!(samp)
    # Use JuMP to solve the linear programming problem
    ## Make a copy of the solver before each call of the warmup so it is thread safe. Because of that force the use of a single thread.
    m = Model(() -> Gurobi.Optimizer(GRB_ENV))
    set_optimizer_attribute(m, "Threads", 1)
    
    @variable(m, samp.prob.lower[i] <= α[i = 1:size(samp.prob.E, 2)] <= samp.prob.upper[i])
    @constraint(m, samp.prob.E * α + samp.prob.r .== 0)
    @constraint(m, samp.prob.L * α .>= 0)

    # Warm up along each individual dimension to get the Min and Max to form a spanning set
    for sense in (MOI.MIN_SENSE, MOI.MAX_SENSE)
        for i in 1:size(samp.prob.E, 2)
            @objective(m, sense, α[i])
            optimize!(m)
            if termination_status(m) == MOI.INFEASIBLE
                samp.feasible = false
                return :Infeasible
            end
            if termination_status(m) == MOI.OPTIMAL
                push!(samp.warmup, value.(α))
            end
        end
    end
    # Remove redundant search directions
    unique!(samp.warmup)
    samp.n_warmup = length(samp.warmup)
    if samp.n_warmup == 0
        error("Insufficient Direction Vectors Found")
    end
    samp.n_samples = samp.n_warmup
    samp.center = mean(samp.warmup)

    # start at the current estimated center
    samp.prev = copy(samp.center)

    samp.feasible = true
    return :Feasible
end

"""
Sample a new feasible point from the point `sampler.prev` in direction `Δ`.
"""
function hr_step(samp, Δin)
    Δ = copy(Δin)
    prob = samp.prob
    #TODO: should add this to sampler type
    feasibility_tol = 1e-7
    scale_neg = -Inf
    scale_pos = Inf
    for i in eachindex(Δ)
        # If any of the directions are near zero we should not sample in that direction as it will
        # not be in the convex space (i.e. that dimension is fixed and will not vary)
        # if isapprox(Δ[i], 0, atol = feasibility_tol)
        #     # A degenerate direction dimensions, skip this dimension
        #     #@show "Hit and Run search direction is length 0"
        #     continue
        # end
        if abs(Δ[i]) < feasibility_tol
            Δ[i] = 0.0
            continue
        end
        # We need to solve:
        # α + Δ * scale = Bound => scale = (Bound - α) / Δ
        #
        # (1 - bounds_tol):
        # make the bounds a little tighter according to the `bounds_tol` so that we have less
        # of a chance of leaving the boundary due to numerical error
        scale_low = ((1 - samp.bounds_tol) * prob.lower[i] - samp.prev[i]) / Δ[i]
        scale_upp = ((1 - samp.bounds_tol) * prob.upper[i] - samp.prev[i]) / Δ[i]

        for sv in (scale_low, scale_upp)
            if sv < 0
                # Find the largest negative scaling
                if sv > scale_neg
                    scale_neg = sv
                end
            else
                # Find the smallest positive scaling
                if sv < scale_pos
                    scale_pos = sv
                end
            end
        end
    end
    # Inequality bounds
    # Inequalitity terms
    Lx = samp.prob.L * samp.prev
    LΔ = samp.prob.L * Δ

    #TODO: would be nice to have bounds_tol fixed in the sampler type so it can vary for each inequality
    for i in eachindex(LΔ)
        if abs(LΔ[i]) < feasibility_tol
            continue
        end
        #ieq_low = ((1 - samp.bounds_tol) * prob.lower[i] - Lx[i]) / LΔ[i]
        #ieq_upp = ((1 - samp.bounds_tol) * prob.upper[i] - Lx[i]) / LΔ[i]
        sv = (samp.bounds_tol - Lx[i]) / LΔ[i]

        if sv < 0
            # Find the largest negative scaling
            if sv > scale_neg
                scale_neg = sv
            end
        else
            # Find the smallest positive scaling
            if sv < scale_pos
                scale_pos = sv
            end
        end
    end
    #@show scale_neg, scale_pos
    scale = rand(samp.RNG, Uniform(scale_neg, scale_pos))

    return samp.prev + scale * Δ
end

function update!(samp::HRSampler, Δ)
    samp.prev = hr_step(samp, Δ)
    samp.n_samples += 1
end

function update!(samp::ACHRSampler, Δ)
    samp.prev = hr_step(samp, Δ)
    # update the center given the new iteration
    samp.center = ((samp.n_samples * samp.center) / (samp.n_samples + 1) + samp.prev / (samp.n_samples + 1))
    samp.n_samples += 1
end

function hr_sample!(samp::HRSampler, n::Int; thinning::Int = 1, burn_in::Int = -1)
    if burn_in == -1 # what is the correct way to do this in a typed manner?
        burn_in = round(Int, n / 2)
    end
    chain = Vector{eltype(samp.prev)}[]

    for i in 1:n
        # mix in the original warmup points to not get stuck
        Δ = samp.Z * rand(size(samp.Z, 2))
        update!(samp, Δ)
        # thin for memory saving
        if i % thinning == 0 && i > burn_in
            push!(chain, copy(samp.prev))
            #push!(chain, copy(samp.center))
        end
    end

    return VectorOfArray(chain)
end

"""
Generate a set of samples.

this is the basic sampling function for all hit-and-run samplers.

* Arguments
n::Int
    the number of samples that are generated at once

* Notes
Performance of this function linearly depends on the number of reactions in your model and the
thinning factor.

ACHR generates samples by choosing new directions from the sampling space's
center and the warmup points. The implementation used here is the same
as in the Matlab Cobra Toolbox [2]_ and uses only the initial warmup points
to generate new directions and not any other previous iterates. This
usually gives better mixing since the startup points are chosen to span
the space in a wide manner. This also makes the generated sampling chain
quasi-markovian since the center converges rapidly.
Memory usage is roughly in the order of (2 * number reactions)^2
due to the required nullspace matrices and warmup points. So large
models easily take up a few GB of RAM.

* References
[1] Direction Choice for Accelerated Convergence in Hit-and-Run Sampling
   David E. Kaufman Robert L. Smith
   Operations Research 199846:1 , 84-95
   https://doi.org/10.1287/opre.46.1.84
[2] https://github.com/opencobra/cobratoolbox
"""
function hr_sample!(samp::ACHRSampler, n::Int; thinning::Int = 1, burn_in::Int = -1)
    if burn_in == -1
        burn_in = round(Int, n / 2)
    end
    chain = Vector{eltype(samp.prev)}[]

    for i in 1:n
        # mix in the original warmup points to not get stuck
        Δ = sample(samp.warmup) - samp.center
        update!(samp, Δ)
        # thin for memory saving
        if i % thinning == 0 && i > burn_in
            push!(chain, copy(samp.prev))
            #push!(chain, copy(samp.center))
        end
    end

    return VectorOfArray(chain)
end
