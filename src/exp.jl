function all_hyperparams()
    gs::Vector{Real} = logspace(0.1, 2, 5)
    ms::Vector{Real} = logspace(0.1, 2, 5)
    ϵs::Vector{Real} = logspace(0.975, 0.999, 5)
    ηs::Vector{Real} = logspace(1e-3, 1e-1, 3)
    Δs::Vector{Real} = logspace(1e-4, 1e-2, 3)
    return product(gs, ms, ϵs, ηs, Δs)
end

function nlopt_grid_search(
    E::Function, d::Integer;
    num_solutions::Integer=10,
    θ_start::Union{Nothing,Vector}=nothing,
    constraints::Union{Nothing,Vector}=nothing,
    seed::Integer=42,
)::Tuple{Solution,HyperParameters,Real}
    Random.seed!(seed)

    best = nothing
    
    for (i, (g, m, ϵ, η, Δ)) in enumerate(all_hyperparams())
        result = nothing
        try
            print("$i\t")
            result = @time rmc(
                E, d, num_solutions;
                g=g, m=m, ϵ=ϵ, η=η, Δ=Δ,
                θ_start=θ_start,
                constraints=constraints,
                isobjective=true,
            )
        catch e
            @warn e g m ϵ η Δ
            continue
        end

        solution = argmin(s -> E(s.answer), result.solutions)
        val = E(solution.answer)
        if isnothing(best)
            # initial case
            best = (solution, result.hyperparams, val)
        elseif val ≈ best[3] && solution.evaluations < best[1].evaluations
            # special case: more or less found the minimum, but want to see if
            # we did it in fewer evalutions of E
            best = (solution, result.hyperparams, val)
        elseif val < best[3]
            # obvious case: found a better minimum
            best = (solution, result.hyperparams, val)
        end
    end

    @assert all(C -> !iscollision(C, best[1].answer), constraints) "constraints broken (bug?)"

    return best
end
