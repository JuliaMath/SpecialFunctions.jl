using SpecialFunctions, Documenter

makedocs(modules=[SpecialFunctions],
         sitename="SpecialFunctions.jl",
         authors="Jeff Bezanson, Stefan Karpinski, Viral B. Shah, et al.",
         format = Documenter.HTML(; assets = String[]),
         pages=["Home" => "index.md",
                "Overview" => "functions_overview.md",
                "Reference" => "functions_list.md"])

deploydocs(repo="github.com/JuliaMath/SpecialFunctions.jl.git")
