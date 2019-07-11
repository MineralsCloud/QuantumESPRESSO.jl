using Documenter, QuantumESPRESSO

makedocs(;
    modules=[QuantumESPRESSO],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/singularitti/QuantumESPRESSO.jl/blob/{commit}{path}#L{line}",
    sitename="QuantumESPRESSO.jl",
    authors="Qi Zhang <singularitti@outlook.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/singularitti/QuantumESPRESSO.jl",
)
