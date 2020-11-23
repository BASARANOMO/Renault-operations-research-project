struct Emballage
    e::Int
    l::Int

    Emballage(; e, l) = new(e, l)
end

function Base.show(io::IO, emballage::Emballage)
    str = "\nEmballage $(emballage.e)"
    str *= "\n   Longueur $(emballage.l)"
    print(io, str)
end

function lire_emballage(row::String, dims::NamedTuple)::Emballage
    row_split = split(row, r"\s+")
    e = parse(Int, row_split[2]) + 1
    l = parse(Int, row_split[4])
    return Emballage(e = e, l = l)
end