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

function check_stock(usine::Usine)::Bool
    # Constraint (2)
    @assert all(0 .<= usine.s)
    return true
end

function check_stock(fournisseur::Fournisseur)::Bool
    # Constraint (4b)
    @assert all(0 .<= fournisseur.s)
    return true
end

function feasibility(instance::Instance)::Bool
    J, U, F, E = instance.J, instance.U, instance.F, instance.E

    for route in instance.routes
        check_length(route)
        check_load(route, instance)
    end

    for usine in instance.usines
        check_stock(usine)
    end
    for fournisseur in instance.fournisseurs
        check_stock(fournisseur)
    end

    return true
end