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

## Extended Solar System (JPL Horizons)

The extended API adds Pluto, the Moon, the Galilean moons and Titan by
querying the JPL Horizons service. **These calls require network access.**

```julia
using AstroIC, Dates

# Sun + 8 planets + Pluto + Moon + Io/Europa/Ganymede/Callisto + Titan (16 bodies)
bodies = solar_system_full(DateTime(2026, 7, 25))

# A single moon, relative to its parent, in SI units
m = moon(DateTime(2026, 7, 25); body = "Europa")

# Raw Horizons state vector: (x, y, z, vx, vy, vz) in km and km/s
sv = horizons_state_vector("301", DateTime(2026, 7, 25); center = "500@399")
```

Results are cached in memory per `(body, center, step, epoch)`. Call
`horizons_state_vector` directly if you need a state vector in a different
reference center. Positions are heliocentric ecliptic J2000 Cartesian
coordinates, consistent with [`solarsystem`](@ref).

!!! tip "Use [AstroPlot.jl](https://github.com/JuliaAstroSim/AstroPlot.jl) to visualize the inner Solar System"