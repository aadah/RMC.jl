import RMC

using Distributions

begin
    d = 1

    m = 2.5
    dist = MixtureModel([Normal(μ, 1) for μ in [-m, m]])
    E(θ) = logpdf(dist, θ[1])

    RMC.mcmc_grid_search(
        E, d,
        islogenergy=true,
        target=dist,
        saveto="$(@__FILE__).csv"
    )
end
