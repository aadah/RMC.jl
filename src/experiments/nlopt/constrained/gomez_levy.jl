import RMC

begin
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

    solution, hyperparams, min_val = RMC.nlopt_grid_search(
        F, 2,
        θ_start=[-0.9, 0.9], # adversarial start point
        constraints=[C],
        seed=42,
    )

    # solution should be ≈ -1.031 @ [0.089, -0.712]
    @info "result" repr(solution.answer) solution.evaluations min_val hyperparams
end
