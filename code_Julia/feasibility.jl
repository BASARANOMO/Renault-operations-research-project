function check_length(route::Route)::Bool
    @assert nb_stops(route) <= 4
    return true
end

function check_load(route::Route, instance::Instance)::Bool
    q = sum(stop.Q for stop in route.stops)
    l_route = sum(q[e] * instance.emballages[e].l for e = 1:instance.E)
    # Constraint (8)
    @assert l_route <= instance.L

    return true
end

function check_stock(usine::Usine, s::Matrix{Int})::Bool
    # Constraint (2)
    @assert all(0 .<= s)
    return true
end

function check_stock(fournisseur::Fournisseur, s::Matrix{Int})::Bool
    # Constraint (4b)
    @assert all(0 .<= s)
    return true
end

function feasibility(solution::Solution, instance::Instance)::Bool
    J, U, F, E = instance.J, instance.U, instance.F, instance.E

    for route in solution.routes
        check_length(route)
        check_load(route, instance)
    end

    su, sf = compute_stocks(solution, instance)

    for usine in instance.usines
        check_stock(usine, su[:, usine.u, :])
    end
    for fournisseur in instance.fournisseurs
        check_stock(fournisseur, sf[:, fournisseur.f, :])
    end

    return true
end