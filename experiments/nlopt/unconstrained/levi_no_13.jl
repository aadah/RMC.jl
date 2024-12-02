# Source: https://www.nber.org/system/files/working_papers/w26340/w26340.pdf

import RMC

begin
    d = 10

    F(θ) = sin(3π * θ[1])^2 +
           (θ[end] - 1)^2 * (1 + sin(2π * θ[end])^2) +
           sum(
               idx -> (θ[idx] - 1)^2 * (1 + sin(3π * θ[idx+1])^2),
               1:d-1
           ) +
           1

    solution, hyperparams, min_val = RMC.nlopt_grid_search(
        F, d,
        θ_start=RMC.sample_hypersphere_surface( # adversarial start point
            d,
            radius=10,
            center=ones(d)
        ),
        seed=42,
    )

    # solution should be ~= 1 @ ones(d)
    @info "result" repr(solution.answer) solution.evaluations min_val hyperparams
end
