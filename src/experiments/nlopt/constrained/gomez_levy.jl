begin
    # Random.seed!(42)

    num_sol = 10

    F(θ) = 4 * θ[1]^2 -
           2.1 * θ[1]^4 +
           (1 / 3) * θ[1]^6 +
           θ[1] * θ[2] -
           4 * θ[2]^2 +
           4 * θ[2]^4

    # constraints are regulary spaced "pillars"
    C(θ) = 1.5 +
           sin(4π * θ[1]) -
           2 * sin(2π * θ[2])^2

    result = @time rmc(
        F, 2, num_sol,
        eta=1e-5,
        isobjective=true,
        constraints=[C],
        θ_start=[-0.9, 0.9], # adversarial start point
        save_trajectory=true,
    )

    log_summary(result)

    solution = argmin(F, result.solutions)

    # solution should be ~= -1.031 @ [0.089, -0.712]
    @info "solution" repr(solution) F(solution) C(solution)

    plot_objective(F, [C], xmin=-1, xmax=1, ymin=-1, ymax=1)
end
