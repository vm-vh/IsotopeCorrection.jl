using NaturalIsotopeCorrection
using Documenter

DocMeta.setdocmeta!(NaturalIsotopeCorrection, :DocTestSetup, :(using NaturalIsotopeCorrection); recursive=true)

makedocs(;
    modules=[NaturalIsotopeCorrection],
    authors="Vincent M. von Häfen",
    sitename="NaturalIsotopeCorrection.jl",
    format=Documenter.HTML(;
        canonical="https://vm-vh.github.io/NaturalIsotopeCorrection.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/vm-vh/NaturalIsotopeCorrection.jl",
    devbranch="main",
)
