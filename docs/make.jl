using ApproximateGPs
using Documenter

DocMeta.setdocmeta!(ApproximateGPs, :DocTestSetup, :(using ApproximateGPs); recursive=true)

makedocs(;
    modules=[ApproximateGPs],
    authors="Theo Galy-Fajou <theo.galyfajou@gmail.com> and contributors",
    repo="https://github.com/theogf/ApproximateGPs.jl/blob/{commit}{path}#{line}",
    sitename="ApproximateGPs.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://theogf.github.io/ApproximateGPs.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/theogf/ApproximateGPs.jl",
)
