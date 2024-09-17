function all_hyperparams()
    gs::Vector{Real} = logspace(0.1, 2, 5)
    ms::Vector{Real} = logspace(0.1, 2, 5)
    ϵs::Vector{Real} = logspace(0.975, 0.999, 5)
    ηs::Vector{Real} = logspace(1e-3, 1e-1, 3)
    Δs::Vector{Real} = logspace(1e-4, 1e-2, 3)
    return product(gs, ms, ϵs, ηs, Δs)
end

function nlopt_grid_search(
    F::Function, d::Integer;
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
            stats = @timed rmc(
                F, d, num_solutions;
                g=g, m=m, ϵ=ϵ, η=η, Δ=Δ,
                θ_start=θ_start,
                constraints=constraints,
                isobjective=true,
            )
            result = stats.value
            @info "Experiment $i" elapsed_seconds=stats.time
        catch error
            @error "Experiment $i" error g m ϵ η Δ
            continue
        end

        solution = argmin(s -> F(s.answer), result.solutions)
        val = F(solution.answer)
        if isnothing(best)
            # initial case
            best = (solution, result.hyperparams, val)
        elseif val ≈ best[3] && solution.evaluations < best[1].evaluations
            # special case: more or less found the minimum, but want to see if
            # we did it in fewer evalutions of F
            best = (solution, result.hyperparams, val)
        elseif val < best[3]
            # obvious case: found a better minimum
            best = (solution, result.hyperparams, val)
        end
    end

    if !isnothing(constraints)
        @assert all(C -> !iscollision(C, best[1].answer), constraints) "constraints broken (bug?)"
    end

    return best
end
