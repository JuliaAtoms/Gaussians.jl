using Gaussians
using Documenter

makedocs(;
    modules=[Gaussians],
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com> and contributors",
    repo="https://github.com/JuliaAtoms/Gaussians.jl/blob/{commit}{path}#L{line}",
    sitename="Gaussians.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliaatoms.org/Gaussians.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaAtoms/Gaussians.jl",
)
