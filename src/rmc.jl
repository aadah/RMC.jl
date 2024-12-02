θ(q::Vector) = q[1:length(q)-1] # theta
h(q::Vector) = q[end] # height

U(q::Vector, m::Real, g::Real) = m * g * h(q) # potential energy
K(p::Vector, m::Real) = dot(p, p) / (2 * m) # kinetic energy
H(q::Vector, p::Vector, m::Real, g::Real) = U(q, m, g) + K(p, m) # hamiltonian (total energy)

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

struct Solution
    answer::Vector
    evaluations::Integer # number of evaluations of the target function before finding solution
end

mutable struct Parabola
    q::Vector # start position
    p::Vector # start momentum
    t::Union{Real,Missing} # time of travel
end

struct Result
    hyperparams::HyperParameters # hyperparameters
    accepted::Vector # accepted samples
    rejected::Vector # rejected samples
    trajectory::Vector{Parabola} # complete time-indexed trajectory
    solutions::Vector{Solution} # potential answers to non-linear optimization problem
    evaluations::Integer # number of evaluations of the target function for whole simulation
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
    constraints::Union{Nothing,Vector},
    parabola::Parabola,
    g::Real,
    m::Real,
    Δ::Real;
    tol::Real=1e-6
)
    if isnothing(constraints)
        boundaries = [surface]
    else
        boundaries = [surface, constraints...]
    end

    # double time elapsed until we find a position past the point of collision
    t_start, t_end = 0, Δ
    q_end = position(t_end, parabola, g, m)
    while isempty(findall(B -> iscollision(B, q_end), boundaries))
        t_start, t_end = t_end, 2 * t_end
        q_end = position(t_end, parabola, g, m)
        @assert all(!isinf, q_end) "numerical overflow"
    end

    # binary search in time to find the true point of collision
    while t_end - t_start > tol
        t_mid = t_start + (t_end - t_start) / 2 # NOTE: more stable than (a+b)/2
        q_mid = position(t_mid, parabola, g, m)
        if any(B -> iscollision(B, q_mid), boundaries)
            t_end = t_mid
        else
            t_start = t_mid
        end
    end

    # NOTE: the assumption at this point is that the tolerance is small enough
    # that there is only one boundary (surface or constraint) that is being
    # violated. if there is actually more than one, it's like hitting a "corner"
    # where the two or more boundaries meet. since the tolerance is so small, we
    # practically don't care and arbitrarily pick the first one
    idx = findfirst(B -> iscollision(B, position(t_end, parabola, g, m)), boundaries)
    B = boundaries[idx]

    # take t_start (time right before collision) as the solution time
    return t_start, B
end

function rmc(
    E::Function,  # distribution or function of interest
    d::Integer,   # dimension of parameter space
    n::Integer;   # number of desired solutions (or samples)
    g::Real=1.0,  # gravity
    m::Real=1.0,  # mass
    ϵ::Real=0.99, # coefficient of restitution
    η::Real=1e-2, # refresh threshold
    Δ::Real=1e-3, # initial step size for temporal search
    θ_start::Union{Nothing,Vector}=nothing,
    constraints::Union{Nothing,Vector}=nothing,
    isobjective::Bool=false,
    islogenergy::Bool=false, # ignored if isobjective==true
    save_trajectory::Bool=false,
    count_by_samples::Bool=false,
    max_evaluations::Integer=1000000,
)::Result
    accepted, rejected = [], []
    trajectory, solutions = [], []
    evaluations = 0

    # TODO: Enforce invariants and log warnings.
    hyperparams = HyperParameters(g, m, ϵ, η, Δ)

    s = if isobjective
        # `isobjective == true` if E is simply a function to minimize.
        E
    else
        if islogenergy
            # `islogenergy == true` if E is already the log, so just negate.
            θ -> θ |> E |> -
        else
            # default: S(θ) = -log(E(θ))
            θ -> θ |> E |> log |> -
        end
    end

    # Wrap the call of our joint/objective/etc. so that we count every time
    # it is evaluated.
    function S(θ)
        evaluations += 1
        return s(θ)
    end

    # Define our surface constraint. This ordering of the difference is such
    # that a collision is when it is zero or negative.
    surface(q) = h(q) - S(θ(q))

    # Prepare our boundary constraints (if any).
    constraints = isnothing(constraints) ? [] : map(Constraint, constraints)

    # Initialization.
    θ_start = isnothing(θ_start) ? randn(d) : θ_start
    q, p = refresh_qp(θ_start, m, S)

    while (count_by_samples ? length(accepted) : length(solutions)) < n
        @assert evaluations <= max_evaluations "exceeded allowed number of function evaluations ($max_evaluations)"

        # Catch numerical issues and end early.
        if any(isnan, q) || any(isnan, p) || any(isinf, q) || any(isinf, p)
            @warn "early stop due to numerical instability" q p
            break
        end

        # Create hyperparabolic segment with q and p at t=0. `missing` means we
        # have yet to find the time it collides with the surface or constraint.
        parabola = Parabola(q, p, missing)

        # Find the collision time and boundary, which could either be the
        # surface or a constraint.
        parabola.t, hyperplane = temporalsearch(surface, constraints, parabola, g, m, Δ)

        if save_trajectory
            push!(trajectory, parabola)
        end

        # Get the position and momentum at time of collision.
        q = position(parabola, g, m)
        p = momentum(parabola, g, m)

        θ_i = θ(q) # Our candidate sample.

        p_before = p
        p = reflect(p_before, hyperplane, q) # Bounce off of surface or constraint.

        # Catch special case where the bounce is off a constraint.
        if hyperplane isa Constraint
            # We don't remove energy from the system, nor consider this location
            # as a possible solution, so we simple continue.
            continue
        # Flip a biased coin to determine whether to accept the candidate.
        elseif accept_bounce(p_before, p)
            push!(accepted, θ_i)
        # Track rejection if non-constraint bounce not accepted.
        else
            push!(rejected, θ_i)
        end

        p *= ϵ # Simulate entropic loss of kinetic energy.

        # If the particle has "puttered out", then we say the trajectory has run
        # its course, and we refresh to begin a new one from which to take
        # samples. Whether the sample was accepted or not from a sampling
        # perspective, the sample is a potential solution from the NL opt view.
        if K(p, m) < η
            push!(solutions, Solution(θ_i, evaluations))
            q, p = refresh_qp(θ_i, m, S)
        end
    end

    result = Result(
        hyperparams, # Hyperparameters for the simulation.
        accepted,    # The accepted samples.
        rejected,    # The rejected samples.
        trajectory,  # All the hyperparabolic trajectories followed.
        solutions,   # Possible solution points that minimized our target.
        evaluations  # Total number of evaluations of the target function.
    )

    return result
end
