function all_hyperparams(mcmc::Bool=false)
    gs::Vector{Real} = logspace(0.1, 2, 5)
    ms::Vector{Real} = logspace(0.1, 2, 5)
    ϵs::Vector{Real} = logspace(0.975, 0.999, 5)
    ηs::Vector{Real} = logspace(1e-3, 1e-1, 3)
    Δs::Vector{Real} = logspace(1e-4, 1e-2, 3)
    if mcmc
        push!(ϵs, 1) # don't remove energy from system
        push!(ηs, Inf) # refresh after every single bounce
    end
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
            @info "Experiment $i" elapsed_seconds = stats.time
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

function mcmc_grid_search(
    E::Function, d::Integer;
    num_samples::Integer=1000,
    num_chains::Integer=5, # recommendation: at least 4
    islogenergy::Bool=false,
    random_start_fn::Function=randn,
    seed::Integer=42,
    target::Union{Nothing,UnivariateDistribution}=nothing,
    metric::Function=r -> -r[:rhat_median],
    saveto::Union{Nothing,String}=nothing,
    pval_cutoff::Real=0.05,
)::Dict{Symbol,Any}
    Random.seed!(seed)

    best = nothing

    df = DataFrame()

    for (i, (g, m, ϵ, η, Δ)) in enumerate(all_hyperparams(true))
        chains = Matrix[]
        stats = @timed for j in 1:num_chains
            result = nothing
            try
                result = rmc(
                    E, d, num_samples;
                    g=g, m=m, ϵ=ϵ, η=η, Δ=Δ,
                    islogenergy=islogenergy,
                    count_by_samples=true,
                    θ_start=random_start_fn(d),
                )
                samples = result.accepted
                push!(chains, collect(hcat(samples...)'))
            catch error
                @error "Experiment $i [chain $j]" error g m ϵ η Δ
            end
        end

        if isempty(chains)
            @error "Experiment $i" error = "no chains" g m ϵ η Δ
            continue
        end

        row = diagnostics(chains)
        row[:num_chains] = length(chains)
        row[:run] = i
        row[:gravity] = g
        row[:mass] = m
        row[:restitution] = ϵ
        row[:eta] = η
        row[:step_size] = Δ

        if !isnothing(target)
            # add pval for test statistic (univariate only)
            chains = [chain[:, 1] for chain in chains]
            row[:pval] = ks.(chains, target)
            row[:pval_median] = median(row[:pval])
        end

        push!(df, NamedTuple(row), cols=:union)

        @info "Experiment $i" metric(row) elapsed_seconds = stats.time

        if isnothing(best)
            # initial case
            best = row
        elseif row[:num_chains] == num_chains && metric(row) > metric(best) && get(row, :pval_median, 1.0) >= pval_cutoff
            # conditions to be met:
            # - we were able to run all chains w/o error
            # - the test stat pval doesn't reject our null hypothesis
            # - finally, the metric of interest is better
            best = row
            @info "Current best" best
        end
    end

    if !isnothing(saveto)
        CSV.write(saveto, df)
    end

    return best
end
