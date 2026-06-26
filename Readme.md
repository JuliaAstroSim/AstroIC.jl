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

or directly from git:

```julia
using Pkg; Pkg.add("https://github.com/JuliaAstroSim/AstroIC.jl")
```

To test the package:

```julia
]test AstroIC
```

## Documentation

- [**Dev**][docs-dev-url] &mdash; *documentation of the in-development version.*

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://juliaastrosims.github.io/AstroIC.jl/dev

For beginners, it is highly recommended to read the [documentation of PhysicalParticles.jl](https://juliaastrosim.github.io/PhysicalParticles.jl/dev/).

## What's included

### Initial condition configurations

| Config type         | Model                                                    |
| ------------------- | -------------------------------------------------------- |
| `PlummerStarCluster`| Plummer (1911) spherical cluster                         |
| `ExponentialDisc`   | Exponential disk with optional central hole              |
| `Bulge`             | Sersic-like bulge with axis ratio `q` and power-law `α`  |
| `GasCloud`          | Cartesian grid SPH gas cloud                             |
| `SphericalSystem`   | Generic spherical system from a discrete/analytic mass shell |

### Particle manipulation tools

`setpos` / `setvel` shift a particle system to a target position / velocity,
`addpos` / `addvel` translate by a fixed offset.

### Real-world data loaders

Pre-bundled CSV / DAT tables from the literature:

- Milky Way satellites ([Battaglia+ 2022](https://www.aanda.org/articles/aa/full_html/2022/08/aa43935-22))
- Ultra-faint dwarf galaxies
- SPARC late- and early-type galaxy rotation curves
- Milky Way rotation curves (Eilers 2019, Mróz 2019, Wang 2021)
- Massive dwarf rotation curves (Cooke 2022)

## Usage

### Plummer star cluster

```julia
using AstroIC
using PhysicalParticles, UnitfulAstro

# Keyword-only constructor; defaults to 1000 particles, 0.01 kpc virial radius,
# 1e5 Msun, Newtonian gravity, uAstro units.
config = PlummerStarCluster(
    collection   = STAR,
    NumSamples   = 100,
    VirialRadius = 0.010u"kpc",
    TotalMass    = 1.0e5u"Msun",
    model        = AstroSimBase.Newton(),
)

# MaxRadius restricts the sampling region (warns if smaller than VirialRadius).
particles = generate(config; MaxRadius = 0.1u"kpc")

# Pass `uSI` to get SI units instead of the default uAstro.
particles_si = generate(config, uSI)
```

### Exponential disk (with optional central hole)

```julia
config = ExponentialDisc(
    collection   = STAR,
    NumSamples   = 10_000,
    TotalMass    = 1.0e10u"Msun",
    ScaleRadius  = 2.0u"kpc",
    ScaleHeight  = 0.02u"kpc",
    HoleRadius   = 0.0u"kpc",   # set > 0 for a central hole
)

particles = generate(config)
```

Pass a `RotationCurve = (radii, velocities)` pair to assign realistic
rotational velocities instead of zero.

### Bulge

```julia
config = Bulge(
    collection  = STAR,
    NumSamples  = 5_000,
    TotalMass   = 8.57u"Msun",
    ScaleRadius = 0.075u"kpc",
    CutRadius   = 2.1u"kpc",
    q           = 0.5,
    α           = 1.8,
)

particles = generate(config)
```

### Spherical system from a mass shell

```julia
using DataFrames

# Either an analytic shell-mass function:
shell_mass_func = x -> 1.0e8u"Msun/kpc^3" / (x / (2.0u"kpc")) /
                              (1 + x / (2.0u"kpc"))^2 * 4π * x^2

config = SphericalSystem(STAR, 1000, shell_mass_func)
particles = generate(config; MaxRadius = 6.0u"kpc")
```

Or a discrete `(r, m)` table as a `Dict` or `DataFrame` with columns `:r` /
`:m`.

### Gas cloud (SPH)

```julia
config = GasCloud(
    collection   = GAS,
    Radius       = 20u"kpc",
    rho0         = 1250u"Msun/kpc^3",
    T            = 300u"K",
    ParticleMass = Constant().m_p,
    Nx           = 11, Ny = 11, Nz = 11,
)

particles = generate(config)  # default uAstro
```

### Solar system

```julia
using Dates
using AstroIC

particles = solarsystem(now())  # SI units; returns 8 planets + the Sun
```

### Setting system center & velocity

```julia
data = generate(PlummerStarCluster(NumSamples = 100))

# Shift the system so its median position equals the given point.
setpos(data, PVector(100.0, 100.0, 100.0, u"kpc"))

# Set mass-weighted mean velocity.
setvel(data, PVector(100.0, 100.0, 100.0, u"kpc/Gyr"))
```

### Loading real-world data

```julia
# Milky Way satellites (Battaglia+ 2022): 61 entries
df = load_data_MW_satellites()

# SPARC late-type galaxy rotation curves (Lelli+ 2016)
df = load_SPARC_LTGs_data()
rc = load_SPARC_LTGs_RC()

# Milky Way rotation curves from three independent measurements
df_eilers = load_MW_RC_Eilers2019()
df_mroz   = load_MW_RC_Mroz2019()
df_w21    = load_MW_RC_DS_W21()
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
