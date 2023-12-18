function sample_hypersphere_surface(d; radius=1, center=nothing)
    v = rand(Dirichlet(d, 1))
    v = sqrt.(v)
    v *= radius
    v = isnothing(center) ? v : v + center
    return v
end
