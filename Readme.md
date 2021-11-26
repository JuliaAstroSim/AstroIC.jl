# AstroIC.jl

Initial condition generator for astrophysical simulations

[![codecov](https://codecov.io/gh/JuliaAstroSim/AstroIC.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaAstroSim/AstroIC.jl)
[![][docs-dev-img]][docs-dev-url]

## Installation

```julia
]add AstroIC
```

or

```julia
using Pkg; Pkg.add("AstroIC")
```

or

```julia
using Pkg; Pkg.add("https://github.com/JuliaAstroSim/AstroIC.jl")
```

To test the Package:
```julia
]test AstroIC
```

## Documentation

- [**Dev**][docs-dev-url] &mdash; *documentation of the in-development version.*

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://juliaastrosim.github.io/AstroIC.jl/dev

For beginners, it is highly recommended to read the [documentation of PhysicalParticles.jl](https://juliaastrosim.github.io/PhysicalParticles.jl/dev/).

## Usage

### Generating a star cluster using Plummer model

```julia
using AstroIC
using PhysicalParticles, UnitfulAstro

## First define a config. Keywords are necessary since the config type is immutable
config = PlummerStarCluster(
    collection = STAR,
    NumSamples = 100,
    VirialRadius = 0.010u"kpc",
    TotalMass = 1.0e5u"Msun",
    model = AstroIC.Newton(),
)

## Now generate particles. MaxRadius restricts the sampling region.
particles = generate(config, MaxRadius = 0.1u"kpc")

# Default units is uAstro, to use SI units:
particles = generate(config, uSI)
```

### Generating a gas cloud

```julia
using AstroIC
using PhysicalParticles, UnitfulAstro

config = GasCloud(
    collection = GAS,
    Radius = 20u"kpc",
    rho0 = 1250u"Msun/kpc^3",
    T = 300u"K",
    ParticleMass = Constant().m_p,
    Nx = 11,
    Ny = 11,
    Nz = 11,
)

particles = generate(config)
```

### Generating our solar system

```julia
using Dates
using AstroIC
particles = solarsystem(now())  # SI units
```

## Package ecosystem

- Basic data structure: [PhysicalParticles.jl](https://github.com/JuliaAstroSim/PhysicalParticles.jl)
- File I/O: [AstroIO.jl](https://github.com/JuliaAstroSim/AstroIO.jl)
- Initial Condition: [AstroIC.jl](https://github.com/JuliaAstroSim/AstroIC.jl)
- Parallelism: [ParallelOperations.jl](https://github.com/JuliaAstroSim/ParallelOperations.jl)
- Trees: [PhysicalTrees.jl](https://github.com/JuliaAstroSim/PhysicalTrees.jl)
- Meshes: [PhysicalMeshes.jl](https://github.com/JuliaAstroSim/PhysicalMeshes.jl)
- Plotting: [AstroPlot.jl](https://github.com/JuliaAstroSim/AstroPlot.jl)
- Simulation: [ISLENT](https://github.com/JuliaAstroSim/ISLENT)