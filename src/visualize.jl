function plot_surface_2d(S; origin=(-1, -1), dx=2, dy=2, fig_size=(6.9, 4))
    # w, h = fig_size
    # set_default_plot_size(w * inch, h * inch)

    Gadfly.with_theme(:dark) do
        return plot(
            z=(x, y) -> S([x, y]),
            xmin=[origin[1]], xmax=[origin[1] + dx],
            ymin=[origin[2]], ymax=[origin[2] + dy],
            Geom.contour,
        )
    end
end

function plot_samples_1d(result::Result, S; from=-1, to=1, fig_size=(6.9, 4))
    # w, h = fig_size
    # set_default_plot_size(w * inch, h * inch)

    samples = [samp[1] for samp in result.accepted]
    Gadfly.with_theme(:dark) do
        p = plot(θ -> S([θ]), from, to, color=[colorant"gold"])
        push!(p, layer(
            x=samples,
            Geom.density,
            color=[colorant"red"],
            linestyle=[:dot],
        ))
        return p
    end
end

function plot_samples_1d_v2(result::Result; dim=1, transform=identity, fig_size=(6.9, 4))
    # w, h = fig_size
    # set_default_plot_size(w * inch, h * inch)

    samples = [transform(samp[dim]) for samp in result.accepted]
    Gadfly.with_theme(:dark) do
        p = plot(
            x=samples,
            Geom.density,
            color=[colorant"red"],
            linestyle=[:dot],
        )
        return p
    end
end

function plot_trajectory(result::Result, S; from=-1, to=1, fig_size=(6.9, 4))
    # w, h = fig_size
    # set_default_plot_size(w * inch, h * inch)

    all_positions = []
    for parabola in result.trajectory
        append!(
            all_positions,
            parabola(
                result.hyperparams.g,
                result.hyperparams.m,
                result.hyperparams.Δ,
            ))
    end
    T = hcat(all_positions...)
    Gadfly.with_theme(:dark) do
        p = plot(θ -> S([θ]), from, to, color=[colorant"gold"])
        push!(p, layer(
            x=T[1, :], y=T[2, :],
            Geom.path,
            linestyle=[:dot],
            color=[colorant"green"],
            Theme(line_width=0.1mm),
        ))
        x = [θ[1] for θ in result.accepted]
        y = [S(θ) for θ in result.accepted]
        push!(p, layer(
            x=x,
            y=y,
            Geom.point,
            color=[colorant"red"],
            Theme(point_shapes=[Shape.diamond]),
        ))
        return p
    end
end

function plot_objective(
    F::Function,
    Cs::Union{Nothing,Vector}=nothing;
    xmin=-1, xmax=1,
    ymin=-1, ymax=1,
    w=6
)
    dx, dy = xmax - xmin, ymax - ymin
    ratio = dy / dx
    h = w * ratio
    set_default_plot_size(w * inch, h * inch)

    f(x, y) = F([x, y])

    Cs = isnothing(Cs) ? [] : Cs
    function a(x, y)
        if any(C -> C([x, y]) <= 0, Cs)
            @info "ZERO"
            return 0
        end
        @info "ONE"
        return 1
    end

    return plot(
        z=f,
        # layer(z=a, Geom.contour,
        #     xmin=[xmin], xmax=[xmax],
        #     ymin=[ymin], ymax=[ymax],),
        alpha=a,
        # z=f.(
        #     collect(xmin:Δ:xmax),
        #     collect(ymin:Δ:ymax)'
        # ),
        # alpha=a.(
        #     collect(xmin:Δ:xmax),
        #     collect(ymin:Δ:ymax)'
        # ),
        xmin=[xmin], xmax=[xmax],
        ymin=[ymin], ymax=[ymax],
        Guide.xlabel("θ[1]"),
        Guide.ylabel("θ[2]"),
        Guide.xticks(ticks=nothing),
        Guide.yticks(ticks=nothing),
        Geom.contour,
    )
end
