# Bulge

The bulge density profile is

```math
\rho(R, z) \propto \frac{\exp\!\big(-r'^2 / r_\mathrm{cut}^2\big)}{(1 + r' / r_s)^\alpha},
\quad
r' = \sqrt{R^2 + (z / q)^2},
```

with axis ratio `q` (oblate for `q < 1`) and inner slope `α`.

## Example

```@repl bulge
using AstroIC
using PhysicalParticles, UnitfulAstro

config = Bulge(
    collection  = STAR,
    NumSamples  = 5000,
    TotalMass   = 8.57u"Msun",
    ScaleRadius = 0.075u"kpc",
    CutRadius   = 2.1u"kpc",
    q           = 0.5,
    α           = 1.8,
)

particles = generate(config)
length(particles)
```

See [`Bulge`](@ref) and [`generate(::Bulge)`](@ref) in the
[Library reference](@ref Methods) for full signatures and keyword
arguments.

!!! tip "Use [AstroPlot.jl](https://github.com/JuliaAstroSim/AstroPlot.jl) to visualize"