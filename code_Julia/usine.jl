struct Usine
    u::Int
    v::Int
    coor::Tuple{Float64,Float64}

    cs::Vector{Int}
    s0::Vector{Int}

    b⁺::Matrix{Int}
    r::Matrix{Int}

    Usine(; u, v, coor, cs, s0, b⁺, r) = new(u, v, coor, cs, s0, b⁺, r)
end

function Base.show(io::IO, usine::Usine)
    str = "\nUsine $(usine.u)"
    str *= "\n   Sommet $(usine.v)"
    str *= "\n   Coordonnées $(usine.coor)"
    str *= "\n   Coûts stock $(usine.cs)"
    str *= "\n   Stock initial $(usine.s0)"
    str *= "\n   Libération journalière $(usine.b⁺)"
    str *= "\n   Stock maximal journalier $(usine.r)"
    print(io, str)
end

function lire_usine(row::String, dims::NamedTuple)::Usine
    row_split = split(row, r"\s+")
    u = parse(Int, row_split[2]) + 1
    v = parse(Int, row_split[4]) + 1
    coor = (parse(Float64, row_split[7]), parse(Float64, row_split[6]))
    k = 8

    cs = Vector{Int}(undef, dims.E)
    s0 = Vector{Int}(undef, dims.E)

    k += 1
    for e = 1:dims.E
        cs[e] = parse(Int, row_split[k+3])
        s0[e] = parse(Int, row_split[k+5])
        k += 6
    end

    b⁺ = Matrix{Int}(undef, dims.E, dims.J)
    r = Matrix{Int}(undef, dims.E, dims.J)

    k += 1
    for j = 1:dims.J
        k += 2
        for e = 1:dims.E
            b⁺[e, j] = parse(Int, row_split[k+3])
            r[e, j] = parse(Int, row_split[k+5])
            k += 6
        end
    end

    return Usine(u = u, v = v, coor = coor, cs = cs, s0 = s0, b⁺ = b⁺, r = r)
end