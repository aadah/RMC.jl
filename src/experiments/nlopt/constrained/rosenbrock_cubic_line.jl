begin
    # Random.seed!(42)

    num_sol = 10

    F(θ) = 100 * (θ[2] - θ[1]^2)^2 + (1 - θ[1])^2

    # constrained to a cubic and a line
    C1(θ) = 1 +
            θ[2] -
            (θ[1] - 1)^3
    C2(θ) = 2 -
            θ[1] -
            θ[2]

    result = @time rmc(
        F, 2, num_sol,
        eta=1e-5,
        isobjective=true,
        constraints=[C1, C2],
        θ_start=zeros(2),
    )

    log_summary(result)

    solution = argmin(F, result.solutions)

    # solution should be ~= 0 @ [1, 1]
    @info "solution" repr(solution) F(solution) C1(solution) C2(solution)
end
