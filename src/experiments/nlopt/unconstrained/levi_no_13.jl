# Source: https://www.nber.org/system/files/working_papers/w26340/w26340.pdf

begin
    # Random.seed!(42)

    num_sol = 1

    d = 2
    F(θ) = sin(3π * θ[1])^2 +
           (θ[end] - 1)^2 * (1 + sin(2π * θ[end])^2) +
           sum(
               idx -> (θ[idx] - 1)^2 * (1 + sin(3π * θ[idx+1])^2),
               1:d-1
           ) +
           1
    answer = ones(d)

    result = @time rmc(
        F, d, num_sol,
        η=1e-5,
        Δ=1, # this seems to help "teleport" through valleys
        isobjective=true,
        θ_start=sample_hypersphere_surface( # adversarial start point
            d,
            radius=10,
            center=answer
        ),
    )

    log_summary(result)

    solution = argmin(F, result.solutions)

    # solution should be ~= 1 @ ones(d)
    @info "solution" repr(solution) F(solution)

    if d == 2
        plot_objective(F, xmin=-10, xmax=10, ymin=-10, ymax=10)
    end
end
