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

function logspace(start::Real, stop::Real, n::Integer=2)
    if n == 2
        return [start, stop]
    end
    s = log10(start)
    e = log10(stop)
    step = (e - s) / (n - 1)
    mid = exp10.((s+step):step:(e+1e-10))
    return [start, mid[1:end-1]..., stop]
end
