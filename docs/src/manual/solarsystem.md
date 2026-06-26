# Solar System

[`solarsystem`](@ref) returns the eight planets plus the Sun at a given date
in SI units. Positions are computed via [AstroLib.jl](https://github.com/JuliaAstro/AstroLib.jl)'s
`helio` (heliocentric coordinates), velocities are derived from a one-second
finite difference.

See [`solarsystem`](@ref) in the [Library reference](@ref Methods) for the
full signature (accepts either a `DateTime` or a Julian day number).

## Example

```@repl solarsystem
using AstroIC
using Dates

particles = solarsystem(now())
length(particles)         # 9 — eight planets + the Sun
first(particles.Pos)      # Mercury's heliocentric position
```

You can also pass a Julian date directly:

```julia
using AstroIC
using Dates

particles = solarsystem(jdcnv(DateTime(2024, 1, 1)))
```

!!! tip "Use [AstroPlot.jl](https://github.com/JuliaAstroSim/AstroPlot.jl) to visualize the inner Solar System"