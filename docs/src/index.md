# AstroIC.jl

*Generates initial conditions for astrophysical simulations.*

Source code: [https://github.com/JuliaAstroSim/AstroIC.jl](https://github.com/JuliaAstroSim/AstroIC.jl)

## Package feature

Initial-condition generators for the standard components of a galaxy-like
simulation:

- Spherical cluster: [`PlummerStarCluster`](@ref)
- Galactic disk: [`ExponentialDisc`](@ref) (with optional central hole)
- Galactic bulge: [`Bulge`](@ref)
- SPH gas cloud: [`GasCloud`](@ref)
- Generic spherical system from a mass shell: [`SphericalSystem`](@ref)
- Solar system snapshot: [`solarsystem`](@ref)

Particle-system helpers:

- Position & velocity translation: [`addpos`](@ref), [`addvel`](@ref)
- Centering & frame setting: [`setpos`](@ref), [`setvel`](@ref)
- Cartesian grid generation: [`gridpoints`](@ref)
- Thermal speed: [`vmean`](@ref)

Real-world observational data loaders (pre-bundled with the package):

- Milky Way satellite galaxies
- Ultra-faint dwarf galaxies
- SPARC late- and early-type galaxy rotation curves
- Milky Way rotation curves (Eilers 2019, Mróz 2019, Wang 2021)
- Massive dwarf rotation curves (Cooke 2022)

## Quick start

```@repl index
using AstroIC
using PhysicalParticles, UnitfulAstro

config = PlummerStarCluster(
    NumSamples   = 1000,
    VirialRadius = 0.010u"kpc",
    TotalMass    = 1.0e5u"Msun",
)

particles = generate(config; MaxRadius = 0.050u"kpc")
length(particles)
```

## Manual outline

```@contents
Pages = [
    "manual/guide.md",
    "manual/plummer.md",
    "manual/disk.md",
    "manual/bulge.md",
    "manual/spherical.md",
    "manual/gascloud.md",
    "manual/solarsystem.md",
    "manual/tools.md",
    "manual/data.md",
]
Depth = 2
```

## Library reference

```@contents
Pages = ["lib/Types.md", "lib/Methods.md"]
Depth = 2
```

## Package ecosystem

- Basic data structure: [PhysicalParticles.jl](https://github.com/JuliaAstroSim/PhysicalParticles.jl)
- File I/O: [AstroIO.jl](https://github.com/JuliaAstroSim/AstroIO.jl)
- Initial Condition: [AstroIC.jl](https://github.com/JuliaAstroSim/AstroIC.jl)
- Parallelism: [ParallelOperations.jl](https://github.com/JuliaAstroSim/ParallelOperations.jl)
- Trees: [PhysicalTrees.jl](https://github.com/JuliaAstroSim/PhysicalTrees.jl)
- Meshes: [PhysicalMeshes.jl](https://github.com/JuliaAstroSim/PhysicalMeshes.jl)
- Plotting: [AstroPlot.jl](https://github.com/JuliaAstroSim/AstroPlot.jl)
- [AstroNbodySim.jl](https://github.com/JuliaAstroSim/AstroNbodySim.jl) — gravitational N-body simulations, glue layer, and parallel runtime
- Simulation of wave dark matter: [WaveDM.jl](https://github.com/JuliaAstroSim/WaveDM.jl)
