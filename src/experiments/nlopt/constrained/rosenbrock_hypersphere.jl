begin
    # Random.seed!(42)

    num_sol = 10

    d = 10
    F(θ) = sum(
        idx -> 100 * (θ[idx+1] - θ[idx]^2)^2 + (1 - θ[idx])^2,
        1:d-1
    )

    # constrained to a hypersphere of radius √d
    C(θ) = d - sum(idx -> θ[idx]^2, 1:d)

    result = @time rmc(
        F, d, num_sol,
        η=1e-5,
        isobjective=true,
        constraints=[C],
        θ_start=zeros(d),
    )

    log_summary(result)

    solution = argmin(F, result.solutions)

    # solution should be ~= 0 @ ones(d)
    @info "solution" repr(solution) F(solution) C(solution)
end
