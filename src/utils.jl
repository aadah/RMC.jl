function sample_hypersphere_surface(d; radius=1, center=nothing)
    v = rand(Dirichlet(d, 1))
    v = sqrt.(v)
    v *= radius
    v = isnothing(center) ? v : v + center
    v .*= rand([1, -1], d)
    return v
end

function log_summary(result::Result)
    num_evaluations = result.evaluations
    mean_sample = string(mean(result.accepted))
    num_accepted = length(result.accepted)
    num_rejected = length(result.rejected)
    num_solutions = length(result.solutions)
    mean_samples_per_trajectory = num_accepted / num_solutions
    reject_rate = "$(round(100 * (num_rejected / (num_accepted+num_rejected)), digits=2))%"
    @info "RMC summary" num_evaluations mean_sample num_accepted num_rejected mean_samples_per_trajectory reject_rate num_solutions
end

logspace(start, stop, step=1) = exp10.(start:step:stop)
