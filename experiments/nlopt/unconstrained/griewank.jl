# Source: https://www.nber.org/system/files/working_papers/w26340/w26340.pdf

import RMC

begin
    d = 10

    F(θ) = 2 +
           sum(
        idx -> θ[idx]^2 / 200,
        1:d
    ) -
           prod(
        idx -> cos(θ[idx] / sqrt(idx)),
        1:d
    )

    solution, hyperparams, min_val = RMC.nlopt_grid_search(
        F, d,
        θ_start=RMC.sample_hypersphere_surface( # adversarial start point
            d,
            radius=50,
            center=zeros(d)
        ),
        seed=42,
    )

    # solution should be ~= 1 @ zeros(d)
    @info "result" repr(solution.answer) solution.evaluations min_val hyperparams
end