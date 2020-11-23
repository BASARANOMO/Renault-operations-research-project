struct Fournisseur
    f::Int
    v::Int
    coor::Tuple{Float64,Float64}

    cs::Vector{Int}
    cexc::Vector{Int}
    s0::Vector{Int}

    b⁻::Matrix{Int}
    r::Matrix{Int}

    Fournisseur(; f, v, coor, cs, cexc, s0, b⁻, r) = new(f, v, coor, cs, cexc, s0, b⁻, r)
end

function Base.show(io::IO, fournisseur::Fournisseur)
    str = "\nFournisseur $(fournisseur.f)"
    str *= "\n   Sommet $(fournisseur.v)"
    str *= "\n   Coordonnées $(fournisseur.coor)"
    str *= "\n   Coûts stock $(fournisseur.cs)"
    str *= "\n   Coûts expédition carton $(fournisseur.cexc)"
    str *= "\n   Stock initial $(fournisseur.s0)"
    str *= "\n   Consommation journalière $(fournisseur.b⁻)"
    str *= "\n   Stock maximal journalier $(fournisseur.r)"
    print(io, str)
end

function lire_fournisseur(row::String, dims::NamedTuple)::Fournisseur
    row_split = split(row, r"\s+")
    f = parse(Int, row_split[2]) + 1
    v = parse(Int, row_split[4]) + 1
    coor = (parse(Float64, row_split[7]), parse(Float64, row_split[6]))
    k = 8

    cs = Vector{Int}(undef, dims.E)
    cexc = Vector{Int}(undef, dims.E)
    s0 = Vector{Int}(undef, dims.E)

    k += 1
    for e = 1:dims.E
        cs[e] = parse(Int, row_split[k+3])
        cexc[e] = parse(Int, row_split[k+5])
        s0[e] = parse(Int, row_split[k+7])
        k += 8
    end

    b⁻ = Matrix{Int}(undef, dims.E, dims.J)
    r = Matrix{Int}(undef, dims.E, dims.J)

    k += 1
    for j = 1:dims.J
        k += 2
        for e = 1:dims.E
            b⁻[e, j] = parse(Int, row_split[k+3])
            r[e, j] = parse(Int, row_split[k+5])
            k += 6
        end
    end

    return Fournisseur(
        f = f,
        v = v,
        coor = coor,
        cs = cs,
        cexc = cexc,
        s0 = s0,
        b⁻ = b⁻,
        r = r,
    )
end