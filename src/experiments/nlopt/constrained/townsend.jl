begin
    # Random.seed!(42)

    num_sol = 10

    F(θ) = -cos((θ[1] - 0.1) * θ[2])^2 -
           θ[1] * sin(3 * θ[1] + θ[2])

    atan2(θ) = θ[1] / (θ[2] + norm(θ, 2))
    function C(θ)
        t = atan2(θ)
        return (2 * cos(t) -
                0.5 * cos(2 * t) -
                0.25 * cos(3 * t) -
                0.125 * cos(4 * t))^2 +
               (2 * sin(t))^2 -
               θ[1]^2 -
               θ[2]^2
    end

    result = @time rmc(
        F, 2, num_sol,
        η=1e-5,
        isobjective=true,
        constraints=[C],
        θ_start=[-1, -1.5], # adversarial start point
        save_trajectory=true,
    )

    log_summary(result)

    solution = argmin(F, result.solutions)

    # solution should be ~= -2.023 @ [2.005, 1.194]
    @info "solution" repr(solution) F(solution) C(solution)
end
