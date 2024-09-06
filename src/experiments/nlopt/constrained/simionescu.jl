begin
    F(θ) = 0.1 * θ[1] * θ[2]

    r_T, r_S, n = 1, 0.2, 8
    C(θ) = (r_T + r_S * cos(n * atan(θ[1] / θ[2])))^2 -
           θ[1]^2 -
           θ[2]^2

    solution, hyperparams, min_val = nlopt_grid_search(
        F, d,
        θ_start=[-0.5, 0.5], # adversarial start point
        constraints=[C],
        seed=42,
    )

    # solution should be ~= -0.072 @ [±0.848, ∓0.848]
    @info "result" repr(solution.answer) solution.evaluations min_val hyperparams
end
