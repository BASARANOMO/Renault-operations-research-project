mutable struct Solution
    R::Int
    routes::Vector{Route}

    Solution(; R, routes) = new(R, routes)
end

function Base.show(io::IO, solution::Solution)
    str = "Solution"
    str *= "\n   Nb de routes: $(solution.R)"
    str *= "\n"
    print(io, str)
    for route in solution.routes
        println(io, route)
    end
end

function lire_solution(path::String)
    sol = open(path) do file
        readlines(file)
    end

    R = parse(Int, split(sol[1], r"\s+")[2])
    routes = [lire_route(sol[1+r]) for r = 1:R]

    return Solution(R = R, routes = routes)
end

function compute_stocks(solution::Solution, instance::Instance)
    J, U, F, E = instance.J, instance.U, instance.F, instance.E

    # Constraints (10.a) and (10.b)

    z⁻ = zeros(Int, E, U, J)
    z⁺ = zeros(Int, E, F, J)

    for e = 1:E, j = 1:J
        for route in solution.routes
            for u = 1:U
                z⁻[e, u, j] += pickup(route, u = u, e = e, j = j)
            end
            for f = 1:F
                z⁺[e, f, j] += delivery(route, f = f, e = e, j = j)
            end
        end
    end

    # Constraints (1) and (5)

    b⁺ = collect(instance.usines[u].b⁺[e, j] for e = 1:E, u = 1:U, j = 1:J)
    b⁻ = collect(instance.fournisseurs[f].b⁻[e, j] for e = 1:E, f = 1:F, j = 1:J)

    su0 = collect(instance.usines[u].s0[e] for e = 1:E, u = 1:U)
    sf0 = collect(instance.fournisseurs[f].s0[e] for e = 1:E, f = 1:F)

    su = Array{Int,3}(undef, E, U, J)
    sf = Array{Int,3}(undef, E, F, J)

    for e = 1:E, j = 1:J
        for u = 1:U
            su[e, u, j] = (j == 1 ? su0[e, u] : su[e, u, j-1]) + b⁺[e, u, j] - z⁻[e, u, j]
        end
        for f = 1:F
            sf[e, f, j] =
                max(0, (j == 1 ? sf0[e, f] : sf[e, f, j-1]) - b⁻[e, f, j]) + z⁺[e, f, j]
        end
    end

    return su, sf
end