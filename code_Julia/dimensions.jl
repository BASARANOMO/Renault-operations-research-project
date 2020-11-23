function lire_dimensions(row::String)::NamedTuple
    row_split = split(row, r"\s+")
    return (
        J = parse(Int, row_split[2]),
        U = parse(Int, row_split[4]),
        F = parse(Int, row_split[6]),
        E = parse(Int, row_split[8]),
        L = parse(Int, row_split[10]),
        γ = parse(Int, row_split[12]),
        ccam = parse(Int, row_split[14]),
        cstop = parse(Int, row_split[16]),
    )
end