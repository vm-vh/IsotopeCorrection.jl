using Documenter, Literate, IsotopeCorrection

examples = sort(filter(x -> endswith(x, ".jl"), readdir(joinpath(@__DIR__, "src"), join = true)))

for example in examples
    Literate.markdown(
        example,
        joinpath(@__DIR__, "src"),
        repo_root_url = "https://github.com/vm-vh/IsotopeCorrection.jl/blob/main",
    )
end

example_mds = first.(splitext.(basename.(examples))) .* ".md"

withenv("COLUMNS" => 150) do
    makedocs(
        modules = [IsotopeCorrection],
        clean = false,
        format = Documenter.HTML(
            ansicolor = true,
            canonical="https://vm-vh.github.io/IsotopeCorrection.jl/dev/",
        ),
        sitename = "IsotopeCorrection.jl",
        linkcheck = false,
        pages = [["Background" => "background.md", "Reference" => "index.md"];example_mds],
    )
end

deploydocs(
    repo = "github.com/vm-vh/IsotopeCorrection.jl.git",
    target = "build",
    branch = "gh-pages",
    push_preview = false,
)
