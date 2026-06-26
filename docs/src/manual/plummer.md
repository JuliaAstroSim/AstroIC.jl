# Plummer

The Plummer (1911) model is a classic analytical equilibrium for a spherical
star cluster. AstroIC supports both **Newtonian** and **MOND** (Milgrom 1983)
gravity models.

```math
\rho(r) = \frac{3 M}{4 \pi r_s^3} \Big(1 + \frac{r^2}{r_s^2}\Big)^{-5/2}
```

## Example

```@repl plummer
using AstroIC
using PhysicalParticles, UnitfulAstro

# config
config = PlummerStarCluster(
    collection   = STAR,
    NumSamples   = 1000,
    VirialRadius = 0.010u"kpc",
    TotalMass    = 1.0e5u"Msun",
    model        = AstroSimBase.Newton(),
)

# generate
particles = generate(config; MaxRadius = 0.050u"kpc")
length(particles)
```

See [`PlummerStarCluster`](@ref) and [`generate(::PlummerStarCluster)`](@ref)
in the [Library reference](@ref Methods) for full signatures.

## Keyword reference

| Keyword        | Default       | Description |
| -------------- | ------------- | ----------- |
| `collection`   | `STAR`        | Particle collection type |
| `NumSamples`   | `1000`        | Number of particles |
| `VirialRadius` | `0.010u"kpc"` | Plummer scale radius |
| `TotalMass`    | `1.0e5u"Msun"`| Total mass of the cluster |
| `model`        | `Newton()`    | `AstroSimBase.Newton` or `MOND1983Milgrom` |

!!! tip "Use [AstroPlot.jl](https://github.com/JuliaAstroSim/AstroPlot.jl) to visualize"

## Switching to MOND

```julia
config = PlummerStarCluster(;
    NumSamples   = 1000,
    VirialRadius = 0.010u"kpc",
    TotalMass    = 1.0e5u"Msun",
    model        = MOND1983Milgrom(),
)
particles = generate(config)
```

`MaxRadius` defaults to ``5 \times \text{VirialRadius}``; particles sampled
outside that sphere are rejected and re-sampled.