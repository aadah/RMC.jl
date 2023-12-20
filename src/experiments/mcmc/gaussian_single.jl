begin
    # Random.seed!(42)

    set_default_plot_size(10 * inch, 8 * inch)

    d = 1
    μ = -10
    σ = 3
    dist = MvNormal(fill(μ, d), I(d) * σ)
    P(θ) = pdf(dist, θ)
    S(θ) = -log(P(θ))
    from, to = μ - 15, μ + 15

    num_sol = 50

    result = @time rmc(
        P, d, num_sol,
        g=6e-3,
        m=1.0,
        ϵ=0.995,
        η=1e-5,
        Δ=1e-3,
    )

    log_summary(result)
    @info "std" std(result.accepted)[1]
    
    vstack(
        # plot_trajectory(result, S, from=from, to=to),
        plot_samples_1d(result, P, from=from, to=to)
        # plot(x=1:length(result.speeds), y=result.speeds, Geom.line)
    )
end
