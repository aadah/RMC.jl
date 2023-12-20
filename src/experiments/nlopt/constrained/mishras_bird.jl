begin
    # Random.seed!(42)

    num_sol = 10

    F(θ) = (θ[1] - θ[2])^2 +
           cos(θ[1]) * exp((1 - sin(θ[2]))^2) +
           sin(θ[2]) * exp((1 - cos(θ[1]))^2)

    # constrained to radius 5 circle centered at [-5,-5]
    C(θ) = 25 -
           (θ[1] + 5)^2 -
           (θ[2] + 5)^2

    result = @time rmc(
        F, 2, num_sol,
        g=1e-2,
        η=1e-5,
        isobjective=true,
        constraints=[C],
        θ_start=[-8, -8], # adversarial start point
        save_trajectory=true,
    )

    log_summary(result)

    solution = argmin(F, result.solutions)

    # solution should be ~= -106.76 @ [-3.13, -1.58]
    @info "solution" repr(solution) F(solution) C(solution)

    plot_objective(F, [C], xmin=-10, xmax=0, ymin=-10, ymax=0)
end
