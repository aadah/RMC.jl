# Exp 2: mixture of varied gaussians

begin
	Random.seed!(42)
	set_default_plot_size(10 * inch, 8 * inch)

	d = 1
	dists = [MvNormal(fill(p*10+randn()*p, d), I(d) * (1/p)) for p in 1:3]
	P(θ) = sum([pdf(dist, θ) for dist in dists]) / 3
	S(θ) = -log(P(θ))
	from, to = 0, 40

    # mysampler = getsamples
    mysampler = getsamples_v2
    result = @time mysampler(
		P, d, 100,
		g=1e-2,
		m=1,
		ϵ=0.99,
		η=1e-5,
		Δ=1e-1,
        # keep_all_bounces=true,
        count_trajectories=true,
	)
	
	log_summary(result)

	vstack(
		# plot_trajectory(result, S, from=from, to=to),
		plot_samples_1d(result, P, from=from, to=to)
		# plot(x=1:length(result.speeds), y=result.speeds, Geom.line)
	)
end
