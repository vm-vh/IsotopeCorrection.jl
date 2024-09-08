#push!(LOAD_PATH,"../src/")
using NaturalIsotopeCorrection
using Documenter

#DocMeta.setdocmeta!(NaturalIsotopeCorrection, :DocTestSetup, :(using NaturalIsotopeCorrection); recursive=true)

makedocs(;
    modules=[NaturalIsotopeCorrection],
    clean = false,
    authors="Vincent M. von Häfen",
    sitename="NaturalIsotopeCorrection.jl",
    format = Documenter.HTML(
        ansicolor = true,
        canonical="https://vm-vh.github.io/NaturalIsotopeCorrection.jl/stable/",
    ),
    linkcheck = false,
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(
    repo = "github.com/vm-vh/NaturalIsotopeCorrection.jl.git",
    target = "build",
    branch = "gh-pages",
    push_preview = false,
)
