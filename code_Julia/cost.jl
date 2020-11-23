function truck_cost(route::Route, instance::Instance)::Int
    return instance.ccam
end

function stops_cost(route::Route, instance::Instance)::Int
    return instance.cstop * nb_stops(route)
end

function km_cost(route::Route, instance::Instance)::Int
    return instance.γ * nb_km(route, instance)
end

function cost_single_truck(route::Route, instance::Instance)::Int
    return truck_cost(route, instance) +
           stops_cost(route, instance) +
           km_cost(route, instance)
end

function cost(route::Route, instance::Instance)::Int
    return route.x * cost_single_truck(route, instance)
end

function stock_cost(usine::Usine, s::Matrix{Int}; j::Int)::Int
    E, J = size(s)
    return sum(usine.cs[e] * max(0.0, s[e, j] - usine.r[e, j]) for e = 1:E)
end

function cost(usine::Usine, s::Matrix{Int}; j::Int)::Int
    return stock_cost(usine, s, j = j)
end

function cost(usine::Usine, s::Matrix{Int})::Int
    E, J = size(s)
    return sum(cost(usine, s, j = j) for j = 1:J)
end

function stock_cost(fournisseur::Fournisseur, s::Matrix{Int}; j::Int)::Int
    E, J = size(s)
    c = sum(fournisseur.cs[e] * max(0.0, s[e, j] - fournisseur.r[e, j]) for e = 1:E)
    return c
end

function expedition_cost(fournisseur::Fournisseur, s::Matrix{Int}; j::Int)::Int
    E, J = size(s)
    c = sum(
        fournisseur.cexc[e] *
        max(0.0, fournisseur.b⁻[e, j] - (j == 1 ? fournisseur.s0[e] : s[e, j-1]))
        for e = 1:E
    )
    return c
end

function cost(fournisseur::Fournisseur, s::Matrix{Int}; j::Int)::Int
    return stock_cost(fournisseur, s, j = j) + expedition_cost(fournisseur, s, j = j)
end

function cost(fournisseur::Fournisseur, s::Matrix{Int})::Int
    E, J = size(s)
    return sum(cost(fournisseur, s, j = j) for j = 1:J)
end

function cost_verbose(solution::Solution, instance::Instance)::Int
    U, F, J = instance.U, instance.F, instance.J
    usines, fournisseurs = instance.usines, instance.fournisseurs

    su, sf = compute_stocks(solution, instance)

    c = 0

    for u = 1:U
        println("Usine $u")
        for j = 1:J
            cujs = stock_cost(usines[u], su[:, u, :], j = j)
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
            cfjs = stock_cost(fournisseurs[f], sf[:, f, :], j = j)
            cfje = expedition_cost(fournisseurs[f], sf[:, f, :], j = j)
            c += cfjs + cfje
            println("   Jour $j")
            println("      Coût stock: $cfjs")
            println("      Coût expédition: $cfje")
            println("      Coût total: " * string(cfjs + cfje))
        end
    end
    cf = c - cu

    println()

    for route in solution.routes
        println("Route $(route.r), jour $(route.j)")
        crt = truck_cost(route, instance)
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

function cost(solution::Solution, instance::Instance; verbose::Bool = false)::Int
    if verbose
        return cost_verbose(solution, instance)
    end

    U, F, J = instance.U, instance.F, instance.J
    usines, fournisseurs = instance.usines, instance.fournisseurs

    su, sf = compute_stocks(solution, instance)

    c = 0.0

    for u = 1:U
        c += cost(usines[u], su[:, u, :])
    end

    for f = 1:F
        c += cost(fournisseurs[f], sf[:, f, :])
    end

    for route in solution.routes
        c += cost(route, instance)
    end

    return c

end