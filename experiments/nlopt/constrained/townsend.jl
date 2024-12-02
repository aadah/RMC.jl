import RMC

begin
    F(θ) = -cos((θ[1] - 0.1) * θ[2])^2 -
           θ[1] * sin(3 * θ[1] + θ[2])

    function C(θ)
        t = atan(θ[1], θ[2])
        return (2 * cos(t) -
                0.5 * cos(2 * t) -
                0.25 * cos(3 * t) -
                0.125 * cos(4 * t))^2 +
               (2 * sin(t))^2 -
               θ[1]^2 -
               θ[2]^2
    end

    solution, hyperparams, min_val = RMC.nlopt_grid_search(
        F, 2,
        θ_start=[-1, -1.5], # adversarial start point
        constraints=[C],
        seed=42,
    )

    # solution should be ~= -2.023 @ [2.005, 1.194]
    @info "result" repr(solution.answer) solution.evaluations min_val hyperparams
end
