
function nb_stops(route::Route)::Int
    return route.F
end

function nb_km(route::Route, instance::Instance)::Int
    usines, fournisseurs = instance.usines, instance.fournisseurs
    path = [
        usines[route.u].v
        [fournisseurs[stop.f].v for stop in route.stops]
    ]
    return sum(instance.graphe.d[s, t] for (s, t) in zip(path[1:end-1], path[2:end]))
end

function truck_cost(instance::Instance)::Int
    return instance.ccam
end

function stops_cost(route::Route, instance::Instance)::Int
    return instance.cstop * nb_stops(route)
end

function km_cost(route::Route, instance::Instance)::Int
    return instance.γ * nb_km(route, instance)
end

function cost_single_truck(route::Route, instance::Instance)::Int
    return truck_cost(instance) +
           stops_cost(route, instance) +
           km_cost(route, instance)
end

function cost(route::Route, instance::Instance)::Int
    return route.x * cost_single_truck(route, instance)
end

function stock_cost(usine::Usine; j::Int)::Int
    E, J = size(usine.s)
    return sum(usine.cs[e] * max(0.0, usine.s[e, j] - usine.r[e, j]) for e = 1:E)
end

function cost(usine::Usine; j::Int)::Int
    return stock_cost(usine, j = j)
end

function cost(usine::Usine)::Int
    E, J = size(usine.s)
    return sum(cost(usine, j = j) for j = 1:J)
end

function stock_cost(fournisseur::Fournisseur; j::Int)::Int
    E, J = size(fournisseur.s)
    c = sum(
        fournisseur.cs[e] * max(0.0, fournisseur.s[e, j] - fournisseur.r[e, j]) for e = 1:E
    )
    return c
end

function expedition_cost(fournisseur::Fournisseur; j::Int)::Int
    E, J = size(fournisseur.s)
    c = sum(
        fournisseur.cexc[e] * max(
            0.0,
            fournisseur.b⁻[e, j] - (j == 1 ? fournisseur.s0[e] : fournisseur.s[e, j-1]),
        ) for e = 1:E
    )
    return c
end

function cost(fournisseur::Fournisseur; j::Int)::Int
    return stock_cost(fournisseur, j = j) + expedition_cost(fournisseur, j = j)
end

function cost(fournisseur::Fournisseur)::Int
    E, J = size(fournisseur.s)
    return sum(cost(fournisseur, j = j) for j = 1:J)
end

function cost_verbose(instance::Instance)::Int
    U, F, J = instance.U, instance.F, instance.J
    usines, fournisseurs, routes = instance.usines, instance.fournisseurs, instance.routes

    c = 0

    for u = 1:U
        println("Usine $u")
        for j = 1:J
            cujs = stock_cost(usines[u], j = j)
            c += cujs
            println("   Jour $j")
            println("      Coût stock: $cujs")
        end
    end
    cu = c

    println()

    for f = 1:F
        println("Fournisseur $f")
        for j = 1:J
            cfjs = stock_cost(fournisseurs[f], j = j)
            cfje = expedition_cost(fournisseurs[f], j = j)
            c += cfjs + cfje
            println("   Jour $j")
            println("      Coût stock: $cfjs")
            println("      Coût expédition: $cfje")
            println("      Coût total: " * string(cfjs + cfje))
        end
    end
    cf = c - cu

    println()

    for route in routes
        println("Route $(route.r), jour $(route.j)")
        crt = truck_cost(instance)
        crs = stops_cost(route, instance)
        crk = km_cost(route, instance)
        rx = route.x
        c += route.x * (crt + crs + crk)
        println("   Coût camion: $crt")
        println("   Coût arrêts: $crs")
        println("   Coût kilométrique: $crk")
        println("   Nb camions: $rx")
        println("   Coût total: " * string(route.x * (crt + crs + crk)))
    end
    cr = c - cf - cu

    println()
    println("Coût total usines: $cu")
    println("Coût total fournisseurs: $cf")
    println("Coût total routes: $cr")
    println()
    println("Coût total: $c")

    return c
end

function cost(instance::Instance; verbose::Bool = false)::Int
    if verbose
        return cost_verbose(instance)
    end

    c = 0.0
    for usine in instance.usines
        c += cost(usine)
    end
    for fournisseur in instance.fournisseurs
        c += cost(fournisseur)
    end
    for route in instance.routes
        c += cost(route, instance)
    end

    return c

end