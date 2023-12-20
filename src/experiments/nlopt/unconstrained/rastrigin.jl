# Source: https://www.nber.org/system/files/working_papers/w26340/w26340.pdf

begin
    # Random.seed!(42)

    num_sol = 1

    d = 2
    F(θ) = 1 +
           10 * d +
           sum(
               idx -> θ[idx]^2 - 10 * cos(2π * θ[idx]),
               1:d
           )
    answer = zeros(d)

    result = @time rmc(
        F, d, num_sol,
        # g=1e-1,
        η=1e-5,
        Δ=1,
        isobjective=true,
        θ_start=sample_hypersphere_surface( # adversarial start point
            d,
            radius=10,
            center=answer
        ),
    )

    log_summary(result)

    solution = argmin(F, result.solutions)

    # solution should be ~= 1 @ zeros(d)
    @info "solution" repr(solution) F(solution)

    if d == 2
        plot_objective(F, xmin=-5, xmax=5, ymin=-5, ymax=5)
    end
end
