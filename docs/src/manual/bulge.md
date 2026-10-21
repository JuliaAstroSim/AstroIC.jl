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
    TotalMass   = 8.57e9u"Msun",
    ScaleRadius = 0.075u"kpc",
    CutRadius   = 2.1u"kpc",
    q           = 0.5,
    α           = 1.8,
)

particles = generate(config)
length(particles)
```

## With a rotation curve

Pass `RotationCurve = (radii, velocities)` to assign rotational velocities:

```julia
rc_r = [0.5, 1.0, 2.0, 4.0, 6.0] .* u"kpc"
rc_v = [50.0, 100.0, 150.0, 200.0, 220.0] .* u"km/s"
particles = generate(config; RotationCurve = (rc_r, rc_v))
```

## With velocity dispersion

`VelocityDispersion` accepts a constant or a radial function and adds
isotropic Gaussian velocities:

```julia
particles = generate(config; VelocityDispersion = 100.0u"km/s")
particles = generate(config; VelocityDispersion = r -> 100.0u"km/s" * (r / 1.0u"kpc"))
```

See [`Bulge`](@ref) and [`generate(::Bulge)`](@ref) in the
[Library reference](@ref Methods) for full signatures and keyword
arguments.

!!! tip "Use [AstroPlot.jl](https://github.com/JuliaAstroSim/AstroPlot.jl) to visualize"