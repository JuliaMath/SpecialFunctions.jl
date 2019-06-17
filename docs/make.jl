using SpecialFunctions, Documenter

makedocs(
    modules     = [SpecialFunctions],
    sitename    = "SpecialFunctions.jl",
    authors     = "Jeff Bezanson, Stefan Karpinski, Viral B. Shah, et al.",
    pages       = [
        "Home"      => "index.md",
        "Functions" => "special.md",
    ],
)

deploydocs(
    repo    = "github.com/JuliaMath/SpecialFunctions.jl.git",
    target  = "build",
)