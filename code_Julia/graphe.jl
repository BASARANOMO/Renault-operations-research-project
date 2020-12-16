using SparseArrays
using LightGraphs
using ProgressMeter

struct Graphe
    G::SimpleDiGraph
    d::Matrix{Int}

    Graphe(; G, d) = new(G, d)
end

function Base.show(io::IO, graphe::Graphe)
    n, m = nv(graphe.G), ne(graphe.G)
    str = "\nGraphe pondéré avec $n sommets et $m arcs"
    k = min(nv(graphe.G), 10)
    str *= "\n   Distances ($k x $k): " * string(graphe.d[1:k, 1:k])
    print(io, str)
end

function lire_arc(row::String)::NamedTuple
    row_split = split(row, r"\s+")
    v1 = parse(Int, row_split[2]) + 1
    v2 = parse(Int, row_split[3]) + 1
    d = parse(Int, row_split[5])
    return (v1 = v1, v2 = v2, d = d)
end

function lire_graphe(rows::Vector{String}, dims::NamedTuple)::Graphe
    G = SimpleDiGraph(dims.U + dims.F)
    d = zeros(Int, dims.U + dims.F, dims.U + dims.F)
    @showprogress "Reading graph " for row in rows
        a = lire_arc(row)
        add_edge!(G, a.v1, a.v2)
        d[a.v1, a.v2] = a.d
    end
    return Graphe(G = G, d = d)
end