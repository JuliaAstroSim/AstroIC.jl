# Plummer

```@docs
PlummerStarCluster
generate(::PlummerStarCluster)
```

```@repl plummer
using AstroIC
using UnitfulAstro

# config
config = PlummerStarCluster(
    collection = STAR,
    NumSamples = 1000,
    VirialRadius = 0.010u"kpc",
    TotalMass = 1.0e5u"Msun",
    model = AstroSimBase.Newton(),
)

# generate
particles = generate(
    config,
    MaxRadius = 0.050u"kpc",
)
```

!!! tip "Use [AstroPlot.jl](https://github.com/JuliaAstroSim/AstroPlot.jl) to visualize"
