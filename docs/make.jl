using Documenter, Literate, NaturalIsotopeCorrection

examples = sort(filter(x -> endswith(x, ".jl"), readdir(joinpath(@__DIR__, "src"), join = true)))

for example in examples
    Literate.markdown(
        example,
        joinpath(@__DIR__, "src"),
        repo_root_url = "https://github.com/vm-vh/NaturalIsotopeCorrection.jl/blob/main",
    )
end

example_mds = first.(splitext.(basename.(examples))) .* ".md"

withenv("COLUMNS" => 150) do
    makedocs(
        modules = [NaturalIsotopeCorrection],
        clean = false,
        format = Documenter.HTML(
            ansicolor = true,
            canonical="https://vm-vh.github.io/NaturalIsotopeCorrection.jl/dev/",
        ),
        sitename = "NaturalIsotopeCorrection.jl",
        linkcheck = false,
        pages = [["Background" => "background.md", "Reference" => "index.md"];example_mds],
    )
end

deploydocs(
    repo = "github.com/vm-vh/NaturalIsotopeCorrection.jl.git",
    target = "build",
    branch = "gh-pages",
    push_preview = false,
)
