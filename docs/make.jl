using Documenter, NumCME 

makedocs(sitename="𝕹𝖚𝖒ℂ𝕄𝔼",
        authors="Huy Vo",
        doctest=false,
        clean=false,
        modules=[NumCME],        
        format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),        
        pages = [
            "Home" => "index.md",
            "Examples" => Any[
                "examples/telegraph.md",
                "examples/hog1p.md"
            ],
            "API reference" => Any[
                "api/models.md",
                "api/propensities.md",
                "api/propensity_gradients.md",
                "api/space.md",
                "api/vectors.md",
                "api/cmesolve.md",
                "api/cmesenssolve.md"
            ]
        ]
        )