import RMC

using Distributions

begin
    d = 1

    dist = Normal(1, 2)
    E(θ) = logpdf(dist, θ[1])

    best = RMC.mcmc_grid_search(
        E, d,
        islogenergy=true,
        target=dist,
        saveto="$(@__FILE__).csv"
    )

    @info "result" best
end
