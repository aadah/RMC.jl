# Exp 1: single Gaussian

begin
    Random.seed!(44)
    set_default_plot_size(5 * inch, 8 * inch)

    d = 1
    μ = -10
    σ = 3
    dist = MvNormal(fill(μ, d), I(d) * σ)
    P(θ) = pdf(dist, θ)
    S(θ) = -log(P(θ))
    width = 10
    from, to = μ - width, μ + width

    constraints = [
        θ -> θ[1] - (-20), # should always be > -20
        θ -> (-5) - θ[1], # should always be < -5
    ]

    result = @time rmc(
        P, d, 10,
        g=1e-2,
        m=1e-1,
        ϵ=1,
        η=0.001,
        θ_start=[μ - width + 1],
        # constraints=constraints,
        save_trajectory=true,
        count_by_samples=true,
    )

	log_summary(result)
    @info "var" var(result.accepted)[1]

    vstack(
        plot_trajectory(result, S, from=from, to=to),
        plot_samples_1d(result, P, from=from, to=to)
        # plot(x=1:length(result.speeds), y=result.speeds, Geom.line)
    )
end
