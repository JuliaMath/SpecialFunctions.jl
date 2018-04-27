using SpecialFunctions, Documenter

makedocs(
    modules = [SpecialFunctions],
    clean = false,
    format = :html,
    sitename = "SpecialFunctions.jl",
    authors = "Jeff Bezanson, Stefan Karpinski, Viral B. Shah, et al.",
    pages = [
        "Home" => "index.md",
        "Functions" => "special.md",
    ],
)

deploydocs(
    julia = "nightly",
    repo = "github.com/JuliaMath/SpecialFunctions.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
)
