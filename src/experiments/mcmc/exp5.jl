# Exp 1: single Gaussian

begin
    Random.seed!(44)
    set_default_plot_size(10 * inch, 8 * inch)

	α_true, β_true = 4.3, 1.2
	beta_true = Beta(α_true, β_true)

    d = 1
    P(θ) = pdf(beta_true, θ[1])
    S(θ) = -log(P(θ))
    from, to = 0, 1

    cs = [
        RMCConstraint(1, 0, >),
		RMCConstraint(1, 1, <)
    ]

    result = getsamples(
        P, d, 2,
        g=5e-1,
        m=1,
        ϵ=0.99,
        η=0.01,
        Δ=1e-2,
        θ_start=[0.5],
        constraints=cs,
        save_trajectory=true,
        count_trajectories=true,
    )

    @info "results" mean(result.accepted)[1] var(result.accepted)[1]

    vstack(
        plot_trajectory(result, S, from=from, to=to),
        plot_samples_1d(result, P, from=from, to=to)
        # plot(x=1:length(result.speeds), y=result.speeds, Geom.line)
    )
end
