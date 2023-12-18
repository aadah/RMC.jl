# Exp 4: beta (log of) hyperparameter estimation, i.e. no constraint

begin
    Random.seed!(42)
    set_default_plot_size(10 * inch, 8 * inch)

	α_true, β_true = 4.3, 1.2
	beta_true = Beta(α_true, β_true)
	n = 1000
	x = rand(beta_true, n)

    d = 2
    α0, β0 = 1, 1
    α_prior, β_prior = Exponential(α0), Exponential(β0)
    logprior(θ) = logpdf(α_prior, exp(θ[1])) + logpdf(β_prior, exp(θ[2]))
    loglikelihood(θ) = sum(logpdf.(Beta(exp.(θ)...), x))
	logE(θ) = loglikelihood(θ) + logprior(θ)

    result = getsamples(
        logE, d, 2,
        g=1e-2,
        m=1,
        cor=0.99,
        eta=0.001,
        eps=1e-2,
		islogenergy=true,
		count_trajectories=true,
    )

	rmc_ans = mean([exp.(samp) for samp in result.accepted])
	@info "ans" rmc_ans[1] rmc_ans[2]

	vstack(
		plot_samples_1d_v2(result, dim=1, transform=exp),
		plot_samples_1d_v2(result, dim=2, transform=exp)
	)
end
