begin
    # Random.seed!(42)

    num_sol = 10

    F(θ) = 0.1 * θ[1] * θ[2]

    r_T, r_S, n = 1, 0.2, 8
    C(θ) = (r_T + r_S * cos(n * atan(θ[1] / θ[2])))^2 -
           θ[1]^2 -
           θ[2]^2

    result = @time rmc(
        F, 2, num_sol,
        eta=1e-5,
        isobjective=true,
        constraints=[C],
        θ_start=[-0.5, 0.5], # adversarial start point
        save_trajectory=true,
    )

    log_summary(result)

    solution = argmin(F, result.solutions)

    # solution should be ~= -0.072 @ [±0.848, ∓0.848]
    @info "solution" repr(solution) F(solution) C(solution)
end
