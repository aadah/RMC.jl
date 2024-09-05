θ(q::Vector) = q[1:length(q)-1] # theta
h(q::Vector) = q[end] # height

U(q::Vector, m::Real, g::Real) = m * g * h(q) # potential energy
K(p::Vector, m::Real) = dot(p, p) / (2 * m) # kinetic energy
H(q::Vector, p::Vector, m::Real, g::Real) = U(q, m, g) + K(p, m) # hamiltonian (total energy)

speed(p::Vector, m::Real) = norm(p / m) # particle speed, agnostic of direction

function leapfrog(q::Vector, p::Vector, m::Real, g::Real, Δ::Real)
    p[end] -= 0.5 * Δ * m * g
    q += Δ * (p / m)
    p[end] -= 0.5 * Δ * m * g
    return q, p
end

function iscollision(hyperplane::Function, q::Vector)
    return hyperplane(q) <= 0
end

function unit_normal(hyperplane::Function, q::Vector)
    n = gradient(hyperplane, q)[1]
    n /= norm(n)
    return n
end

function reflect(p::Vector, hyperplane::Function, q::Vector)
    n = unit_normal(hyperplane, q)
    p -= 2 * dot(n, p) * n
    return p
end

function accept_bounce(p_before::Vector, p_after::Vector)
    p_before_lat, p_after_lat = θ(p_before), θ(p_after) # we only care about the lateral momentum
    dots = dot(p_before_lat, p_after_lat)
    norms = norm(p_before_lat) * norm(p_after_lat)
    sim = (dots / norms + 1) / 2 # normalize to range [0,1]
    thresh = rand()
    return sim >= thresh
end

function refresh_qp(θ_i::Vector, m::Real, S::Function)
    d = length(θ_i)

    # raise particle to some random height above the surface
    elevation = S(θ_i)
    elevation += rand(Rayleigh(m))

    q = [θ_i..., elevation]
    p = randn(d + 1) * m # "slingshot" in random direction

    return q, p
end

struct Constraint <: Function
    C::Function
end

function (constraint::Constraint)(q)
    return constraint.C(θ(q))
end

struct HyperParameters
    g::Real
    m::Real
    ϵ::Real
    η::Real
    Δ::Real
end

struct Result
    hyperparams::HyperParameters # hyperparameters
    accepted::Vector # accepted samples
    rejected::Vector # rejected samples
    trajectory::Vector # complete time-indexed trajectory
    solutions::Vector # potential answers to non-linear optimization problem
    evaluations::Integer # number of evaluations of the target function
end

mutable struct Parabola
    q::Vector # start position
    p::Vector # start momentum
    t::Union{Real,Missing} # time of travel
end

function position(t::Real, parabola::Parabola, g::Real, m::Real)
    qt = parabola.q + (parabola.p / m) * t
    qt[end] -= (g / 2) * t^2
    return qt
end

function momentum(t::Real, parabola::Parabola, g::Real, m::Real)
    pt = copy(parabola.p)
    pt[end] -= m * g * t
    return pt
end

position(parabola::Parabola, g::Real, m::Real) = position(parabola.t, parabola, g, m)
momentum(parabola::Parabola, g::Real, m::Real) = momentum(parabola.t, parabola, g, m)

# Interpolation
function (parabola::Parabola)(g::Real, m::Real, Δ::Real)
    return [position(t, parabola, g, m) for t in 0:Δ:parabola.t]
end

# much faster discover from reduced number of evaluations of S, avoids MH prob
# that would happen from leapfrog integration, and (relatedly) finds the exact
# collision point. only tradeoff is the "jumping ahead" risk if doubling lands
# us in another "mountain"
function temporalsearch(
    surface::Function,
    parabola::Parabola,
    g::Real,
    m::Real,
    Δ::Real;
    tol::Real=1e-6
)
    # double time elapsed until we find a position past the point of collision
    t_start, t_end = 0, Δ
    q_end = position(t_end, parabola, g, m)
    while !iscollision(surface, q_end)
        t_start = t_end
        t_end *= 2
        q_end = position(t_end, parabola, g, m)
    end

    # binary search in time to find the true point of collision
    while t_end - t_start > tol
        t_mid = (t_start + t_end) / 2
        q_mid = position(t_mid, parabola, g, m)
        if iscollision(surface, q_mid)
            t_end = t_mid
        else
            t_start = t_mid
        end
    end

    # take t_start (time right before collision) as the solution time
    return t_start
end

function temporalsearch(
    constraint::Constraint,
    parabola::Parabola,
    g::Real,
    m::Real;
    tol::Real=1e-6
)
    # quick check if constraint is actually being violated
    t_start, t_end = 0, parabola.t
    q_end = position(t_end, parabola, g, m)
    if !iscollision(constraint, q_end)
        return Inf
    end

    # binary search in time to find the true point of collision
    while t_end - t_start > tol
        t_mid = (t_start + t_end) / 2
        q_mid = position(t_mid, parabola, g, m)
        if iscollision(constraint, q_mid)
            t_end = t_mid
        else
            t_start = t_mid
        end
    end

    # take t_start (time right before collision) as the solution time
    return t_start
end

function rmc(
    E::Function, d::Integer, n::Integer;
    g::Real=1.0,
    m::Real=1.0,
    ϵ::Real=0.99,
    η::Real=1e-2,
    Δ::Real=1e-3,
    θ_start::Union{Nothing,Vector}=nothing,
    constraints::Union{Nothing,Vector}=nothing,
    isobjective::Bool=false,
    islogenergy::Bool=false, # ignored if isobjective==true
    save_trajectory::Bool=false,
    count_by_samples::Bool=false
)::Result
    accepted, rejected = [], []
    trajectory, solutions = [], []
    evaluations = 0

    # TODO: enforce invariants and log warnings
    hyperparams = HyperParameters(g, m, ϵ, η, Δ)

    constraints = isnothing(constraints) ?
                  [] :
                  map(Constraint, constraints)

    s = if isobjective
        # isobjective == true if E is simply a function to minimize
        E
    else
        if islogenergy
            # islogenergy == true if E is already the log, so just negate
            θ -> θ |> E |> -
        else
            # default: S(θ) = -log(E(θ))
            θ -> θ |> E |> log |> -
        end
    end

    # wrap the call of our joint / objective / etc. so that we count every time
    # it is evaluated
    function S(θ)
        evaluations += 1
        return s(θ)
    end

    # this ordering of the difference is such that a collision is when it is
    # zero or negative, i.e. a violation of our "surface" constraint
    surface(q) = h(q) - S(θ(q))

    # initialization
    θ_start = isnothing(θ_start) ? randn(d) : θ_start
    q, p = refresh_qp(θ_start, m, S)

    while (count_by_samples ? length(accepted) : length(solutions)) < n
        # catch numerical issues, end early
        if any(isnan, q) || any(isnan, p) || any(isinf, q) || any(isinf, p)
            @warn "early stop due to numerical instability" q p
            break
        end

        # create hyperparabolic segment with q and p at t=0. `missing` means we
        # have yet to find the time it collides (with the surface or constraint)
        parabola = Parabola(q, p, missing)

        # find the surface collision point. guaranteed to hit the surface, so
        # set `t` on the parabola
        parabola.t = temporalsearch(surface, parabola, g, m, Δ)
        hyperplane = surface

        # TODO: could this end with the original surface constraint being
        # violated? maybe use constraints2 = [surface, constraints...]?
        violated = findall(C -> iscollision(C, position(parabola, g, m)) < 0, constraints)
        while !isempty(violated)
            for constraint in constraints[violated]
                t_c = temporalsearch(constraint, parabola, g, m)
                if t_c < parabola.t
                    parabola.t = t_c
                    hyperplane = constraint # so that reflect will be off constraint
                end
            end
            violated = findall(C -> iscollision(C, position(parabola, g, m)) < 0, constraints)
        end

        if save_trajectory
            push!(trajectory, parabola)
        end

        q = position(parabola, g, m)
        p = momentum(parabola, g, m)

        θ_i = θ(q) # our candidate sample

        p_before = p
        p = reflect(p_before, hyperplane, q) # bounce off of surface or constraint

        if accept_bounce(p_before, p)
            push!(accepted, θ_i)
        else
            push!(rejected, θ_i)
        end

        p *= ϵ # simulate entropic loss of kinetic energy

        # if the particle has "puttered out", then we say the trajectory has run
        # its course, and we refresh to begin a new one from which to take
        # samples. whether the sample was accepted or not from a sampling
        # perspective, the sample is a potential solution from the NL opt view
        if K(p, m) < η
            push!(solutions, θ_i)
            q, p = refresh_qp(θ_i, m, S)
        end
    end

    result = Result(
        hyperparams,
        accepted,
        rejected,
        trajectory,
        solutions,
        evaluations
    )

    log_summary(result)

    return result
end

function log_summary(result::Result)
    num_evaluations = result.evaluations
    mean_sample = string(mean(result.accepted))
    num_accepted = length(result.accepted)
    num_rejected = length(result.rejected)
    num_refreshes = length(result.solutions)
    mean_samples_per_trajectory = num_accepted / num_refreshes
    reject_rate = "$(round(100 * (num_rejected / (num_accepted+num_rejected)), digits=2))%"
    @info "RMC summary" num_evaluations mean_sample num_accepted num_rejected mean_samples_per_trajectory reject_rate num_refreshes
end
