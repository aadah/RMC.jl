# Exp 2: mixture of varied gaussians

begin
	Random.seed!(42)
	set_default_plot_size(10 * inch, 8 * inch)

	d = 1
	dists = [MvNormal(fill(p^2.5, d), I(d) * p) for p in 1:2]
	P(θ) = sum([pdf(dist, θ) for dist in dists]) / 2
	S(θ) = -log(P(θ))
	from, to = -10, 15

	result = getsamples(
		P, d, 10000,
		g=5e-2,
		m=1,
		ϵ=0.99,
		η=0.005,
		Δ=1e-2,
	)
	
	log_summary(result)

	vstack(
		# plot_trajectory(result, S, from=from, to=to),
		plot_samples_1d(result, P, from=from, to=to)
		# plot(x=1:length(result.speeds), y=result.speeds, Geom.line)
	)
end
