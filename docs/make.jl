using Gaussians
using Documenter

DocMeta.setdocmeta!(Gaussians, :DocTestSetup, :(using Gaussians); recursive=true)
makedocs(;
    modules=[Gaussians],
    authors="Stefanos Carlström <stefanos.carlstrom@gmail.com> and contributors",
    repo="https://github.com/JuliaAtoms/Gaussians.jl/blob/{commit}{path}#L{line}",
    sitename="Gaussians.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliaatoms.org/Gaussians.jl",
        assets=String[],
        mathengine = Documenter.MathJax()
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaAtoms/Gaussians.jl",
)
