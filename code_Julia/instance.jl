struct Instance
    J::Int
    U::Int
    F::Int
    E::Int

    L::Int
    γ::Int
    ccam::Int
    cstop::Int

    emballages::Vector{Emballage}
    usines::Vector{Usine}
    fournisseurs::Vector{Fournisseur}
    graphe::Graphe

    Instance(; J, U, F, E, L, γ, ccam, cstop, emballages, usines, fournisseurs, graphe) =
        new(J, U, F, E, L, γ, ccam, cstop, emballages, usines, fournisseurs, graphe)
end

function Base.show(io::IO, instance::Instance)
    str = "\nInstance"
    str *= "\n   Nombre de jours: $(instance.J)"
    str *= "\n   Nombre d'usines: $(instance.U)"
    str *= "\n   Nombre de fournisseurs: $(instance.F)"
    str *= "\n   Nombre de types d'emballages: $(instance.E)"
    print(io, str)
end

function lire_instance(path::String)::Instance
    data = open(path) do file
        readlines(file)
    end

    dims = lire_dimensions(data[1])
    emballages = [lire_emballage(data[1+e], dims) for e = 1:dims.E]
    usines = [lire_usine(data[1+dims.E+u], dims) for u = 1:dims.U]
    fournisseurs = [lire_fournisseur(data[1+dims.E+dims.U+f], dims) for f = 1:dims.F]
    graphe = lire_graphe(data[1+dims.E+dims.U+dims.F+1:end], dims)

    return Instance(
        J = dims.J,
        U = dims.U,
        F = dims.F,
        E = dims.E,
        L = dims.L,
        γ = dims.γ,
        ccam = dims.ccam,
        cstop = dims.cstop,
        emballages = emballages,
        usines = usines,
        fournisseurs = fournisseurs,
        graphe = graphe,
    )
end