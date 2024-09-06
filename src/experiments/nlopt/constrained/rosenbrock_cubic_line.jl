begin
    F(θ) = 100 * (θ[2] - θ[1]^2)^2 + (1 - θ[1])^2

    # constrained to a cubic and a line
    C1(θ) = 1 +
            θ[2] -
            (θ[1] - 1)^3
    C2(θ) = 2 -
            θ[1] -
            θ[2]

    solution, hyperparams, min_val = nlopt_grid_search(
        F, 2,
        θ_start=zeros(2),
        constraints=[C1, C2],
        seed=42,
    )

    # solution should be ~= 0 @ [1, 1]
    @info "result" repr(solution.answer) solution.evaluations min_val hyperparams
end
