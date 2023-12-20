# Source: https://www.nber.org/system/files/working_papers/w26340/w26340.pdf

begin
    # Random.seed!(42)

    num_sol = 10

    d = 2
    F(θ) = 2 +
           sum(
        idx -> θ[idx]^2 / 200,
        1:d
    ) -
           prod(
        idx -> cos(θ[idx] / sqrt(idx)),
        1:d
    )
    answer = zeros(d)

    result = @time rmc(
        F, d, num_sol,
        # g=1e-1,
        η=1e-5,
        # Δ=3, # this seems to help "teleport" through valleys
        isobjective=true,
        θ_start=sample_hypersphere_surface( # adversarial start point
            d,
            radius=50,
            center=answer
        ),
    )

    log_summary(result)

    solution = argmin(F, result.solutions)

    # solution should be ~= 1 @ zeros(d)
    @info "solution" repr(solution) F(solution)

    if d == 2
        plot_objective(F, xmin=-50, xmax=50, ymin=-50, ymax=50)
    end
end
