using Documenter, NumCME 

makedocs(sitename="NumCME",
        authors="Huy Vo",
        doctest=false,
        clean=true,
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

deploydocs(
    repo = "github.com/voduchuy/NumCME.jl.git",
    target = "build",
    deps = Deps.pip("pygments", "mkdocs", "python-markdown-math"),
    devbranch = "main" 
    # ...
)        
