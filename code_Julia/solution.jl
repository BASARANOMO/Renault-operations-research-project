function update_stocks!(usine::Usine, routes::Vector{Route})
    E, J = size(usine.s)

    for e = 1:E, j = 1:J
        usine.z⁻[e, j] = 0
        for route in routes
            usine.z⁻[e, j] += pickup(route, usine, e = e, j = j)
        end
    end

    for e = 1:E, j = 1:J
        usine.s[e, j] =
            (j == 1 ? usine.s0[e] : usine.s[e, j-1]) + usine.b⁺[e, j] - usine.z⁻[e, j]
    end
end

function update_stocks!(fournisseur::Fournisseur, routes::Vector{Route})
    E, J = size(fournisseur.s)

    for e = 1:E, j = 1:J
        fournisseur.z⁺[e, j] = 0
        for route in routes
            fournisseur.z⁺[e, j] += delivery(route, fournisseur, e = e, j = j)
        end
    end

    for e = 1:E, j = 1:J
        fournisseur.s[e, j] =
            max(
                0,
                (j == 1 ? fournisseur.s0[e] : fournisseur.s[e, j-1]) - fournisseur.b⁻[e, j],
            ) + fournisseur.z⁺[e, j]
    end
end


function update_stocks!(instance::Instance, routes::Vector{Route})
    for usine in instance.usines
        update_stocks!(usine, routes)
    end
    for fournisseur in instance.fournisseurs
        update_stocks!(fournisseur, routes)
    end
end

function lire_solution(instance::Instance, path::String)::Instance
    sol = open(path) do file
        readlines(file)
    end

    R = parse(Int, split(sol[1], r"\s+")[2])
    routes = [lire_route(sol[1+r]) for r = 1:R]

    solved_instance = deepcopy(instance)
    solved_instance.R = R
    solved_instance.routes = routes
    update_stocks!(solved_instance, routes)
    
    return solved_instance
end