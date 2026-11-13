using Documenter
using Minimap2

DocMeta.setdocmeta!(Minimap2, :DocTestSetup, :(using Minimap2); recursive=true)

makedocs(;
    modules=[Minimap2],
    sitename="Minimap2.jl",
    authors="Alex Laini",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://alexology.github.io/Minimap2.jl",
        edit_link="main",
    ),
    pages=[
        "Home" => "index.md",
        "Validation" => "validation.md",
        "API" => "api.md",
    ],
    checkdocs=:exports,
)

deploydocs(; repo="github.com/alexology/Minimap2.jl", devbranch="main")
