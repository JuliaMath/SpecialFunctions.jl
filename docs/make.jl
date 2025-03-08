using SpecialFunctions, Documenter
using DocumenterCitations

# `using SpecialFunctions` for all doctests
DocMeta.setdocmeta!(SpecialFunctions, :DocTestSetup, :(using SpecialFunctions); recursive=true)

bib = CitationBibliography(
       joinpath(@__DIR__, "src", "refs.bib");
       style = :authoryear,
)

makedocs(modules=[SpecialFunctions],
         sitename="SpecialFunctions.jl",
         authors="Jeff Bezanson, Stefan Karpinski, Viral B. Shah, et al.",
         format = Documenter.HTML(; assets = String[]),
         pages=["Home" => "index.md",
                "Overview" => "functions_overview.md",
                "Reference" => "functions_list.md"],
         plugins=[bib],
         checkdocs=:exports,
        )

deploydocs(repo="github.com/JuliaMath/SpecialFunctions.jl.git")
