function getsamples(
    E::Function, d::Integer, n::Integer;
    g::Real=1.0,
    m::Real=1.0,
    cor::Real=0.9,
    eta::Real=0.1,
    eps::Real=1e-3,
    θ_start::Union{Nothing,Vector}=nothing,
    islogenergy::Bool=false,
    keep_all_bounces::Bool=false,
    save_trajectory::Bool=false,
    count_trajectories::Bool=false
)::RMCResult
    accepted, rejected = [], []
    trajectory, speeds = [], []
    refreshes = 0

    # islogenergy == true if E is already the log, so we just need to negate.
    S = if islogenergy
        θ -> θ |> E |> -
    else
        θ -> θ |> E |> log |> -
    end

    # This ordering of the difference ensures that the height derivative is
    # possible, i.e. the normal to the gradient will always point upwards. This
    # means the return value is the elevation off the ground, and also that a
    # collision is when it is zero or negative.
    hyperplane(q) = h(q) - S(θ(q))

    θ_start = isnothing(θ_start) ? randn(d) : θ_start
    q, p = refresh_qp(θ_start, m, S)

    H_start = H(q, p, m, g)

    while (count_trajectories ? refreshes : length(accepted)) < n
        # catch numerical issues, end early
        if any(isnan, q) || any(isnan, p) || any(isinf, q) || any(isinf, p)
            break
        end

        # before leapfrog, note location and speed
        if save_trajectory
            push!(trajectory, q)
        end
        push!(speeds, speed(p, m))

        # single leapfrog step
        q_prev = q
        q, p = leapfrog(q, p, m, g, eps)

        if iscollision(hyperplane, q)
            # perform partial time-reversed leapfrog step to estimate true point
            # of collision.
            Δh1 = abs(hyperplane(q))
            Δh2 = abs(hyperplane(q_prev))
            eps_reverse = eps * (Δh1 / (Δh1 + Δh2))
            q, p = leapfrog(q, p, m, g, -eps_reverse)

            # our candidate sample
            θ_i = θ(q)

            # our acceptance prob, per HMC
            H_end = H(q, p, m, g)
            accept = rand() <= min(1, exp(H_start - H_end))

            if accept
                p_prev = p
                q[end] = S(θ_i) # set height to be on S so reflection is accurate
                p = reflect(p_prev, hyperplane, q) # bounce of surface

                if keep_all_bounces || accept_bounce(p_prev, p)
                    push!(accepted, θ_i)
                end

                p *= cor # simulate entropic loss kinetic energy

                # if the particle has "puttered out", then we say the trajectory
                # has run its course, and we refresh to begin a new one from
                # which to take samples
                if K(p, m) < eta
                    q, p = refresh_qp(θ_i, m, S)
                    refreshes += 1
                end
            else
                push!(rejected, θ_i)
                q, p = refresh_qp(θ_i, m, S)
            end

            # collisions either a) are accepted and entropically annealed or b)
            # refreshed completely. either way, we must reset H in order to
            # faithfully do the acceptance probability calculation next
            # collision.
            H_start = H(q, p, m, g)
        end
    end

    result = RMCResult(
        accepted,
        rejected,
        trajectory,
        speeds,
        refreshes,
    )
    log_summary(result)
    return result
end
