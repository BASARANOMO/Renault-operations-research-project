mutable struct RouteStop
    f::Int
    Q::Vector{Int}

    RouteStop(; f, Q) = new(f, Q)
end

mutable struct Route
    r::Int
    j::Int
    x::Int
    u::Int

    F::Int
    stops::Vector{RouteStop}

    Route(; r, j, x, u, F, stops) = new(r, j, x, u, F, stops)
end

function Base.show(io::IO, route::Route)
    str = "Route $(route.r)"
    str *= "\n   Jour $(route.j)"
    str *= "\n   Nb de camions $(route.x)"
    str *= "\n   Usine de départ $(route.u)"
    str *= "\n   Nb d'arrêts $(route.F)"
    for (stoprank, stop) in enumerate(route.stops)
        str *= "\n   Stop $stoprank"
        str *= "\n      Fournisseur $(stop.f)"
        str *= "\n      Livraison $(stop.Q)"
    end
    print(io, str)
end

function lire_route(row::String)::Route
    row_split = split(row, r"\s+")
    r = parse(Int, row_split[2]) + 1
    j = parse(Int, row_split[4]) + 1
    x = parse(Int, row_split[6])
    u = parse(Int, row_split[8]) + 1
    F = parse(Int, row_split[10])

    stops = RouteStop[]

    k = 11
    while k <= length(row_split)
        f = parse(Int, row_split[k+1]) + 1
        k += 2

        Q = Int[]
        while (k <= length(row_split)) && (row_split[k] == "e")
            push!(Q, parse(Int, row_split[k+3]))
            k += 4
        end
        push!(stops, RouteStop(f = f, Q = Q))
    end

    return Route(r = r, j = j, x = x, u = u, F = F, stops = stops)
end

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

function pickup(route::Route, ; u::Int, e::Int, j::Int)::Int
    if (route.j == j) && (route.u == u)
        return route.x * sum(stop.Q[e] for stop in route.stops)
    else
        return 0
    end
end

function delivery(route::Route; f::Int, e::Int, j::Int)::Int
    if route.j == j
        d = 0
        for stop in route.stops
            if stop.f == f
                d += route.x * stop.Q[e]
            end
        end
        return d
    else
        return 0
    end
end