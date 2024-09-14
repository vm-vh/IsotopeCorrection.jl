using Documenter, Literate, NaturalIsotopeCorrection

#DocMeta.setdocmeta!(NaturalIsotopeCorrection, :DocTestSetup, :(using NaturalIsotopeCorrection); recursive=true)

#example = sort(filter(x -> endswith(x, ".jl"), readdir(joinpath(@__DIR__, "src"))))[1]

#Literate.markdown(
#    joinpath(@__DIR__, "src", example),
#    repo_root_url = "https://github.com/vm-vh/NaturalIsotopeCorrection.jl/blob/main",
#)

makedocs(;
    modules=[NaturalIsotopeCorrection],
    clean = false,
    authors="Vincent M. von Häfen",
    sitename="NaturalIsotopeCorrection.jl",
    format = Documenter.HTML(
        ansicolor = true,
        canonical="https://vm-vh.github.io/NaturalIsotopeCorrection.jl/dev/",
    ),
    linkcheck = false,
    pages=[
        "Examples" => "examples.md",
        "Background" => "background.md",
        "Referances" => "index.md",
    ],
)

deploydocs(
    repo = "github.com/vm-vh/NaturalIsotopeCorrection.jl.git",
    target = "build",
    branch = "gh-pages",
    push_preview = false,
)
