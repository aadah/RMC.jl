import RMC

begin
    F(θ) = (θ[1] - θ[2])^2 +
           cos(θ[1]) * exp((1 - sin(θ[2]))^2) +
           sin(θ[2]) * exp((1 - cos(θ[1]))^2)

    # constrained to radius 5 circle centered at [-5,-5]
    C(θ) = 25 -
           (θ[1] + 5)^2 -
           (θ[2] + 5)^2

    solution, hyperparams, min_val = RMC.nlopt_grid_search(
        F, 2,
        θ_start=[-8, -8], # adversarial start point
        constraints=[C],
        seed=42,
    )

    # solution should be ≈ -106.76 @ [-3.13, -1.58]
    @info "result" repr(solution.answer) solution.evaluations min_val hyperparams
end
