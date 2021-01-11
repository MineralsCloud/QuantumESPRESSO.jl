using QuantumESPRESSO
using Documenter

makedocs(;
    modules = [QuantumESPRESSO],
    authors = "Qi Zhang <singularitti@outlook.com>",
    repo = "https://github.com/MineralsCloud/QuantumESPRESSO.jl/blob/{commit}{path}#L{line}",
    sitename = "QuantumESPRESSO.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://MineralsCloud.github.io/QuantumESPRESSO.jl",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/MineralsCloud/QuantumESPRESSO.jl")
