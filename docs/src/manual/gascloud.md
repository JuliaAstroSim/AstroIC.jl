# Gas Cloud

A simple SPH gas cloud sampled on a Cartesian grid clipped to a sphere of
radius `R`. Velocities are drawn from a Maxwell–Boltzmann distribution at
temperature `T`.

## Example

```@repl gascloud
using AstroIC
using PhysicalParticles, UnitfulAstro, PhysicalConstants

config = GasCloud(
    collection   = GAS,
    Radius       = 20u"kpc",
    rho0         = 1250u"Msun/kpc^3",
    T            = 300u"K",
    ParticleMass = Constant().m_p,
    Nx           = 11, Ny = 11, Nz = 11,
)

particles = generate(config)
length(particles)
```

The particle mass per cell is set to ``\rho_0 \cdot R^2 / r^2 \cdot L_x L_y
L_z`` so that the central density is recovered.

See [`GasCloud`](@ref), [`generate(::GasCloud)`](@ref) and [`gridpoints`](@ref)
in the [Library reference](@ref Methods) for full signatures.

!!! tip "Use [AstroPlot.jl](https://github.com/JuliaAstroSim/AstroPlot.jl) to visualize"