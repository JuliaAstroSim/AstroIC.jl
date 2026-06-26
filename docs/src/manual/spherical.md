# Spherical System

A generic spherically symmetric initial condition: pass either an analytic
shell-mass function ``m(r)`` or a discrete table of ``(r_i, m_i)``, and the
generator will rejection-sample ``N`` radii and integrate the total mass
numerically.

`mass_shell` accepts three forms:

- `Function r -> m(r)` with units — analytic profile
- `Dict("r" => r, "m" => m)` — discrete table
- `DataFrame` with columns `:r`, `:m`

## Example: analytic NFW-like profile

```@repl spherical
using AstroIC
using PhysicalParticles, UnitfulAstro

r_s = 2.0u"kpc"
MaxRadius = 3 * r_s
shell_mass_func = x -> 1.0e8u"Msun/kpc^3" / (x / r_s) /
                              (1 + x / r_s)^2 * 4π * x^2

config = SphericalSystem(STAR, 1000, shell_mass_func)
particles = generate(config; MaxRadius)
length(particles)
```

## Example: discrete table

```julia
using DataFrames

r = collect(LinRange(0.001u"kpc", MaxRadius, 50))
m = shell_mass_func.(r)

# As a Dict
config = SphericalSystem(STAR, 1000, Dict("r" => r, "m" => m))
particles = generate(config; MaxRadius)

# As a DataFrame
df = DataFrame(r = r, m = m)
config = SphericalSystem(STAR, 1000, df)
particles = generate(config; MaxRadius)
```

See [`SphericalSystem`](@ref) and [`generate(::SphericalSystem)`](@ref) in the
[Library reference](@ref Methods) for full signatures and accepted forms
of `mass_shell`.

!!! tip "Use [AstroPlot.jl](https://github.com/JuliaAstroSim/AstroPlot.jl) to visualize"