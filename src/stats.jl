import MCMCDiagnosticTools
using HypothesisTests
using Distributions
using Statistics

# Kolmogorov-Smirnov
function ks(samples::Vector, target::UnivariateDistribution)
    return pvalue(ExactOneSampleKSTest(samples, target))
end

function diagnostics(chains::Vector{Matrix})
    chains = stack(chains, dims=2)
    tup = MCMCDiagnosticTools.ess_rhat(chains)
    report = Dict{Symbol,Any}(
        :rhat => tup.rhat, # https://arxiv.org/pdf/1903.08008
        :ess => tup.ess, # effective sample size
        :essr => MCMCDiagnosticTools.ess(chains, relative=true), # effective sample size ratio
    )

    for k in collect(keys(report))
        report[Symbol("$(k)_median")] = median(report[k])
    end

    return report
end
