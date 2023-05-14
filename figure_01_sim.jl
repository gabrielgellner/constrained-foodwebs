using LinearAlgebra
using Distributions
using PyPlot
pygui(true)

λ_stability(M) = maximum(real.(eigvals(M)))

# # Model
# We are solving the model:
# ```math
# Ax + r = 0
# ```
# where
# ```math
# A = \left(
# a_{11} & -a_{12} \\
# a_{21} & 0
# \right)
# ```

# # Check that it satisfies the original model
function A_matrix(α, r, x)
    ## Assuming α = [a[1, 1], a[1, 2]
    A = fill(0.0, 2, 2)
    A[1, 1] = -α[1]
    A[1, 2] = -α[2]
    A[2, 1] = abs(r[2]) / x[1]

    return A
end

# # The Null space of the solution
Z(α12, r, x) = -x[2] * α12 / x[1] + abs(r[1]) / x[1]

function is_feasible(α12, r, x)
    α11 = Z(α12, r, x)
    A = A_matrix([α11, α12], r, x)
    ## Check equality constraints
    ceq = isapprox.(A * x + r, 0, atol = 1e-10)
    ## Check inequality constraints
    ieq = -A[1, 2] >= A[2, 1]
    if all(ceq) && ieq
        return true
    else
        return false
    end
end

# # Generate the community matrix from a sample
M_matrix(α12, r, x) = A_matrix([Z(α12, r, x), α12], r, x) .* x


# # Build all the specific feasibility bounds
function feasibility_par(r, x; upper = 10)
    α12_min = max(0.0, abs(r[2]) / x[1])
    α12_max = min(upper, r[1] / x[2])
    if α12_min > α12_max
        return Float64[]
    end
    return range(α12_min, α12_max, length = 10000)
end


let
    # # Setup the growth / death vector
    r = [1.5, -0.5]
    # # Setup the x[1] value
    x = [1.0, NaN]

    λ1s = []
    max_x2 = 3
    for x_2 in range(1 / max_x2, max_x2, length = 100)
        x[2] = x_2
        α12s = feasibility_par(r, x)
        if length(α12s) > 0
            push!(λ1s, [log(x[2] / x[1]), collect(α12s), [λ_stability(M_matrix(α12, r, x)) for α12 in α12s]])
        end
    end

    summary = fill(0.0, 5, length(λ1s))
    for (i, λ1) in enumerate(λ1s)
        summary[1, i] = λ1[1]
        summary[2, i] = mean(λ1[3])
        summary[3, i] = var(λ1[3])
        summary[4, i] = minimum(λ1[3])
        summary[5, i] = maximum(λ1[3])
    end
    fig, ax1 = plt.subplots()
    ax1.set_xlabel("Log Pyramid Slope (\$log(x_2 / x_1)\$)")
    ax1.set_ylabel("Eigenvalue Stability")
    ax1.plot(summary[1, :], summary[2, :], "r")
    ax1.plot(summary[1, :], summary[4, :], "k")
    ax1.plot(summary[1, :], summary[5, :], "k")

    savefig("figs/figure1-eigstab.png", dpi=300)
end
