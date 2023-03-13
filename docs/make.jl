using Isoplot
using Documenter

DocMeta.setdocmeta!(Isoplot, :DocTestSetup, :(using Isoplot); recursive=true)

makedocs(;
    modules=[Isoplot],
    authors="C. Brenhin Keller",
    repo="https://github.com/JuliaGeochronology/Isoplot.jl/blob/{commit}{path}#{line}",
    sitename="Isoplot.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaGeochronology.github.io/Isoplot.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaGeochronology/Isoplot.jl",
    devbranch = "main",
)
