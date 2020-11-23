
function Base.string(stop::RouteStop)::String
    f = stop.f
    str = "f $(f-1)"
    E = length(stop.Q)
    for e = 1:E
        q = stop.Q[e]
        str *= " e $(e-1) q $q"
    end
    return str
end

function Base.string(route::Route)::String
    r = route.r
    j = route.j
    x = route.x
    u = route.u
    F = route.F
    str = "r $(r-1) j $(j-1) x $x u $(u-1) F $F"
    for stop in route.stops
        str *= " " * string(stop)
    end
    return str
end

function Base.string(solution::Solution)::String
    R = solution.R
    str = "R $R"
    for route in solution.routes
        str *= "\n" * string(route)
    end
    return str
end

function write_sol_to_file(solution::Solution, path::String)::Bool
    open(path, "w") do file
        write(file, string(solution))
    end
    return true
end