# Exponential Disk

The disk density profile is

```math
\rho(R, z) \propto \exp\!\Big(-\frac{|z|}{h_z} - \frac{R}{h_R}\Big),
```

optionally with a central hole (suppressed density for ``R < R_\mathrm{hole}``).

## Example

```@repl disk
using AstroIC
using PhysicalParticles, UnitfulAstro

config = ExponentialDisc(
    collection   = STAR,
    NumSamples   = 5000,
    TotalMass    = 1.0e10u"Msun",
    ScaleRadius  = 2.0u"kpc",
    ScaleHeight  = 0.02u"kpc",
    HoleRadius   = 0.0u"kpc",
)

particles = generate(config; MaxRadius = 6.0u"kpc")
length(particles)
```

See [`ExponentialDisc`](@ref) and [`generate(::ExponentialDisc)`](@ref) in the
[Library reference](@ref Methods) for full signatures and keyword arguments.

## With a rotation curve

Pass `RotationCurve = (radii, velocities)` to assign realistic disk velocities:

```julia
rc_r = [0.5, 1.0, 2.0, 4.0, 6.0] .* u"kpc"
rc_v = [50.0, 100.0, 150.0, 200.0, 220.0] .* u"km/s"
particles = generate(config; RotationCurve = (rc_r, rc_v))
```

The interpolation uses `Dierckx.Spline1D` with keyword arguments `k` (spline
degree, default 2) and `bc` (boundary condition, default `"nearest"`). Set
`rotational_ratio` (default 0.9) to mix in a random isotropic component.

!!! tip "Use [AstroPlot.jl](https://github.com/JuliaAstroSim/AstroPlot.jl) to visualize"