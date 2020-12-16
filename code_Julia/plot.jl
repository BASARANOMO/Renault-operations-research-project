using Plots
import GeoJSON
using GeoInterface

using PyCall
using Conda
Conda.add("folium", channel="conda-forge")

py"""
import folium
"""

function plot_sites(instance::Instance, numbers::Bool = false)
    G = instance.graphe.G
    d = instance.graphe.d
    pl = plot()
    scatter!(
        pl,
        [fournisseur.coor[2] for fournisseur in instance.fournisseurs],
        [fournisseur.coor[1] for fournisseur in instance.fournisseurs],
        marker = (:circle, 5, 0.5, "blue"),
        label = "fournisseurs",
    )

    scatter!(
        pl,
        [usine.coor[2] for usine in instance.usines],
        [usine.coor[1] for usine in instance.usines],
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

    all_y = [site.coor[1] for site in [instance.usines; instance.fournisseurs]]
    all_x = [site.coor[2] for site in [instance.usines; instance.fournisseurs]]
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
    for ((y1, x1), (y2, x2)) in zip(path[1:end-1], path[2:end])
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

function plot_routes(instance::Instance; j::Int)::Plots.Plot
    pl = plot_sites(instance)
    for route in instance.routes
        if route.j == j
            add_route_to_plot!(route, instance, pl)
        end
    end
    plot!(title = "Routes de la solution à j=$j")
    return pl
end

function plot_stocks(usine::Usine)::Plots.Plot
    E, J = size(usine.s)
    pl = plot()
    for e = 1:E
        plot!(pl, 1:J, usine.s[e, :], label = "e=$e", lw = 3, markershape = :circle)
    end
    plot!(
        pl,
        title = "Stocks de l'usine $(usine.u)",
        ylabel = "stock s[e, j]",
        xlabel = "jour j",
        xticks = 1:J,
        xlims = (0.5, J + 0.5),
        ylims = (0, maximum(usine.s) + 1),
        fmt = :png,
    )
    return pl
end

function plot_stocks(
    fournisseur::Fournisseur,
)::Plots.Plot
    E, J = size(fournisseur.s)
    pl = plot()
    for e = 1:E
        plot!(
            pl,
            1:J,
            fournisseur.s[e, :],
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
        ylims = (0, maximum(fournisseur.s[:, :]) + 1),
        fmt = :png,
    )
    return pl
end

function build_geojson_usines(instance::Instance)
    usines_points =
        MultiPoint([Point(usine.coor[1], usine.coor[2]) for usine in instance.usines])
    return FeatureCollection([Feature(usines_points)])
end

function build_geojson_fournisseurs(instance::Instance)
    fournisseurs_points = MultiPoint([
        Point(fournisseur.coor[1], fournisseur.coor[2])
        for fournisseur in instance.fournisseurs
    ])
    return FeatureCollection([Feature(fournisseurs_points)])
end

function build_geojson_routes(instance::Instance)
    usines, fournisseurs = instance.usines, instance.fournisseurs
    linestrings = LineString[]
    for route in instance.routes
        path = [collect(usines[route.u].coor)]
        for stop in route.stops
            push!(path, collect(fournisseurs[stop.f].coor))
        end
        ls = LineString(path)
        push!(linestrings, ls)
    end
    mls = MultiLineString([coordinates(ls) for ls in linestrings])
    return FeatureCollection([Feature(mls)])
end

function plot_folium(instance::Instance; include_fournisseurs::Bool=false)
    py"""
    usines_map = folium.FeatureGroup(name="usines")
    fournisseurs_map = folium.FeatureGroup(name="fournisseurs", show=$include_fournisseurs)
    routes_map = folium.FeatureGroup(name="routes")
    """
    
    for usine in instance.usines
        lat, lon = usine.coor[1], usine.coor[2]
        u = usine.u
        py"""
        folium.Marker(
            location=[$lat, $lon],
            icon=folium.Icon(color='red', icon='cog'),
            popup="Usine {}".format($u)
        ).add_to(usines_map)
        """
    end
    
    for fournisseur in instance.fournisseurs
        lat, lon = fournisseur.coor[1], fournisseur.coor[2]
        f = fournisseur.f
        py"""
        folium.Marker(
            location=[$lat, $lon],
            icon=folium.Icon(color='blue', icon='user'),
            popup="Fournisseur {}".format($f)
        ).add_to(fournisseurs_map)
        """
    end
    
    for route in instance.routes
        r, j, x = route.r, route.j, route.x
        path = [instance.usines[route.u].coor]
        for stop in route.stops
            push!(path, instance.fournisseurs[stop.f].coor)
        end
        py"""
        folium.PolyLine(
            locations=$path,
            color="black",
            popup="Route {} - jour {}".format($r, $j)
        ).add_to(routes_map)
        """
    end
    
    py"""
    m = folium.Map(tiles="Stamen WaterColor", location = [48.85, 2.35], zoom_start = 4)
    folium.TileLayer("OpenStreetMap", show=False, overlay=True).add_to(m)
    usines_map.add_to(m)
    fournisseurs_map.add_to(m)
    routes_map.add_to(m)
    folium.LayerControl().add_to(m)
    """
    
    return py"m"
end