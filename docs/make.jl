using SpecialFunctions, Documenter

makedocs(modules=[SpecialFunctions],
         sitename="SpecialFunctions.jl",
         authors="Jeff Bezanson, Stefan Karpinski, Viral B. Shah, et al.",
         pages=["Home" => "index.md",
                "Overview" => "functions_overview.md",
                "List" => "functions_list.md"])

deploydocs(repo="github.com/JuliaMath/SpecialFunctions.jl.git")
