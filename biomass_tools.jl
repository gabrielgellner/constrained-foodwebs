import Distributions: Uniform
import FoodWebs: trophic_levels


function is_pyramid(pnet, x)
    tls = trophic_levels(pnet)
    for i in eachindex(tls)
        for j in findall(tls[i] .< tls)
            # check all higher trophic levels for having higher biomass
            if x[i] < x[j]
                return false
            end
        end
    end
    return true
end

function is_inverted_pyramid(pnet, x)
    tls = trophic_levels(pnet)
    for i in earchindex(tls)
        for j in findall(tls[i] .< tls)
            # check all higher trophic levels for having lower biomass
            if x[i] > x[j]
                return false
            end
        end
    end
    return true
end

function set_β(pnet, δ, min_β)
    @assert min_β > 0
    # δ is the constant increase / decrease between whole trophic levels
    # so δ > 0 means top heavy,
    #    δ < 0 means eltonian

    tl_scale = abs.(trophic_levels(pnet) .- 1)

    if δ >= 0
        start = min_β
    elseif δ < 0
        # this ensures that:
        # 1. all biomass values will be positive
        # 2. the smallest value will be min_β
        start = min_β - maximum(tl_scale) * δ
    end
    # Fuzz the values by about 20%
    return start .+ tl_scale .* (rand(Uniform(0.8, 1.0), length(tl_scale)) .* δ)
end

function set_r!(conweb, pnet; tp_scale = 1.0)
    tol = 1e-15
    
    R_index = basal_species(pnet)
    for i in 1:pnet.S
        if i in R_index
            conweb.r[i] = rand(Uniform(0, 5))
        else
            conweb.r[i] = -rand(Uniform(0, 1))
        end
    end

    # The structure of E depends on the `r` so we need to update E as well
    conweb.cnet[:, :] = pmat_to_cnet(pnet, conweb.r)
    conweb.E[:, :] = block_diagonal(find_blocks(conweb.cnet, conweb.x))
    return
end
