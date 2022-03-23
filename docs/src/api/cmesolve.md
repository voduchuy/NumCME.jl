# Transient solution of the CME

## FSP variant specification
```@docs
AdaptiveFspSparse
```

## State space adaptation policies
```@autodocs
Modules = [NumCME]
Pages = ["spaceadapters.jl",
        "rstepadapters.jl"
        ]
```

## FSP output format
```@autodocs
Modules = [NumCME]
Pages = [        
        "fspoutput.jl"
]
```
## `solve()` method
```@autodocs
Modules = [NumCME]
Pages = [
    "transientcme/sparse/fspsolve.jl"
]
solve
```
