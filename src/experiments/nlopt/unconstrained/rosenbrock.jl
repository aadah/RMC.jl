# Source: https://www.nber.org/system/files/working_papers/w26340/w26340.pdf

begin
    # Random.seed!(42)

    num_sol = 10

    d = 10
    F(θ) = sum(
        idx -> 100 * (θ[idx+1] - θ[idx]^2)^2 + (1 - θ[idx])^2,
        1:d-1
    ) + 1
    answer = ones(d)

    result = @time rmc(
        F, d, num_sol,
        η=1e-5,
        isobjective=true,
        θ_start=sample_hypersphere_surface( # adversarial start point
            d,
            radius=5,
            center=answer
        ),
    )

    log_summary(result)

    solution = argmin(F, result.solutions)

    # solution should be ~= 1 @ ones(d)
    @info "solution" repr(solution) F(solution)

    if d == 2
        plot_objective(F, xmin=-2, xmax=2, ymin=-1, ymax=2)
    end
end
