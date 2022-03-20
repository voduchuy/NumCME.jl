using Documenter, NumCME 

makedocs(sitename="NumCME documentation",
        authors="Huy Vo",
        doctest=false,
        clean=false,
        modules=[NumCME],
        format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),        
        pages = [
            "Home" => "index.md",
            "API reference" => Any[
                "api/models.md",
                "api/propensities.md",
                "api/propensity_gradients.md",
                "api/vectors.md",
                "api/cmesolve.md",
                "api/cmesenssolve.md"
            ]
        ]
        )