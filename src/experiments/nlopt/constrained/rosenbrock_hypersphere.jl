begin
    d = 10
    F(θ) = sum(
        idx -> 100 * (θ[idx+1] - θ[idx]^2)^2 + (1 - θ[idx])^2,
        1:d-1
    )

    # constrained to a hypersphere of radius √d
    C(θ) = d - sum(idx -> θ[idx]^2, 1:d)

    solution, hyperparams, min_val = nlopt_grid_search(
        F, d,
        θ_start=zeros(d),
        constraints=[C],
        seed=42,
    )

    # solution should be ~= 0 @ ones(d)
    @info "result" repr(solution.answer) solution.evaluations min_val hyperparams
end
