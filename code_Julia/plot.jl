using Plots

function plot_sites(instance::Instance, numbers::Bool = false)
    G = instance.graphe.G
    d = instance.graphe.d
    pl = plot()
    scatter!(
        pl,
        [fournisseur.coor[1] for fournisseur in instance.fournisseurs],
        [fournisseur.coor[2] for fournisseur in instance.fournisseurs],
        marker = (:circle, 5, 0.5, "blue"),
        label = "fournisseurs",
    )

    scatter!(
        pl,
        [usine.coor[1] for usine in instance.usines],
        [usine.coor[2] for usine in instance.usines],
        marker = (:rect, 7, 0.7, "red"),
        label = "usines",
    )

    if numbers
        annotate!(
            pl,
            [
                (
                    fournisseur.coor[1],
                    fournisseur.coor[2],
                    text(string(fournisseur.f), :darkblue, :center),
                ) for fournisseur in instance.fournisseurs
            ],
        )
        annotate!(
            pl,
            [
                (usine.coor[1], usine.coor[2], text(string(usine.u), :darkred, :center))
                for usine in instance.usines
            ],
        )
    end

    all_x = [site.coor[1] for site in [instance.usines; instance.fournisseurs]]
    all_y = [site.coor[2] for site in [instance.usines; instance.fournisseurs]]
    plot!(
        pl,
        grid = false,
        xlabel = "longitude",
        ylabel = "latitude",
        xlim = (minimum(all_x) - 1, maximum(all_x) + 1),
        ylim = (minimum(all_y) - 1, maximum(all_y) + 1),
        title = "Sites de l'instance",
        aspect_ratio = :equal,
        fmt = :png,
    )
    return pl
end

function add_route_to_plot!(route::Route, instance::Instance, pl::Plots.Plot)
    usines, fournisseurs = instance.usines, instance.fournisseurs
    path = [
        usines[route.u].coor
        [fournisseurs[stop.f].coor for stop in route.stops]
    ]
    for ((x1, y1), (x2, y2)) in zip(path[1:end-1], path[2:end])
        plot!(
            pl,
            [x1 + 0.1 * (x2 - x1), x2 - 0.1 * (x2 - x1)],
            [y1 + 0.1 * (y2 - y1), y2 - 0.1 * (y2 - y1)],
            linewidth = 2 * route.x,
            linecolor = "black",
            arrow = arrow(:simple, :triangle, 1.0, 1.0),
            label = nothing,
        )
    end
end

function plot_routes(solution::Solution, instance::Instance; j::Int)::Plots.Plot
    pl = plot_sites(instance)
    for route in solution.routes
        if route.j == j
            add_route_to_plot!(route, instance, pl)
        end
    end
    plot!(title = "Routes de la solution à j=$j")
    return pl
end

function plot_stocks(solution::Solution, instance::Instance, usine::Usine)::Plots.Plot
    su, _ = compute_stocks(solution, instance)

    E, U, J = size(su)
    pl = plot()
    for e = 1:E
        plot!(pl, 1:J, su[e, usine.u, :], label = "e=$e", lw = 3, markershape = :circle)
    end
    plot!(
        pl,
        title = "Stocks de l'usine $(usine.u)",
        ylabel = "stock s[e, j]",
        xlabel = "jour j",
        xticks = 1:J,
        xlims = (0.5, J + 0.5),
        ylims = (0, maximum(su[:, usine.u, :]) + 1),
        fmt = :png,
    )
    return pl
end

function plot_stocks(
    solution::Solution,
    instance::Instance,
    fournisseur::Fournisseur,
)::Plots.Plot
    _, sf = compute_stocks(solution, instance)

    E, F, J = size(sf)
    pl = plot()
    for e = 1:E
        plot!(
            pl,
            1:J,
            sf[e, fournisseur.f, :],
            label = "e=$e",
            lw = 3,
            markershape = :circle,
        )
    end
    plot!(
        pl,
        title = "Stocks du fournisseur $(fournisseur.f)",
        ylabel = "stock s[e, j]",
        xlabel = "jour j",
        xticks = 1:J,
        xlims = (0.5, J + 0.5),
        ylims = (0, maximum(sf[:, fournisseur.f, :]) + 1),
        fmt = :png,
    )
    return pl
end