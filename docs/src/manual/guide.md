# Package Guide

## Installation

From the Julia REPL, type `]` to enter the Pkg REPL mode and run

```julia
pkg> add AstroIC
```

or add from git repository:

```julia
pkg> add https://github.com/JuliaAstroSim/AstroIC.jl
```

Test the package by

```julia
pkg> test AstroIC
```

## Basic usage

The package follows a two-step pattern: first build an immutable **config**
struct describing the system, then call the [`generate`](@ref) dispatcher to
materialize particles.

```@repl guide
using AstroIC
using PhysicalParticles, UnitfulAstro

config = PlummerStarCluster(
    NumSamples   = 1000,
    VirialRadius = 0.010u"kpc",
    TotalMass    = 1.0e5u"Msun",
)
data = generate(config; MaxRadius = 0.050u"kpc")
length(data)
```

By default, particles are returned in `uAstro` units (kpc, kpc/Gyr, Msun). Pass
`uSI` to get SI units instead:

```julia
data_si = generate(config, uSI; MaxRadius = 0.050u"kpc")
```

## Composing a galaxy

The typical workflow is to generate each component separately and concatenate
them into a single particle array:

```@repl guide_compose
using AstroIC
using PhysicalParticles, UnitfulAstro

bulge = generate(Bulge(; NumSamples = 1000, TotalMass = 8.57u"Msun",
                        ScaleRadius = 0.075u"kpc", CutRadius = 2.1u"kpc",
                        q = 0.5, α = 1.8))
disc = generate(ExponentialDisc(; NumSamples = 5000, TotalMass = 1.0e10u"Msun",
                                  ScaleRadius = 2.0u"kpc", ScaleHeight = 0.02u"kpc"))

# Concatenate into a single system
galaxy = vcat(bulge, disc)
length(galaxy)
```

`bulge` and `disc` are `StructArray`s with compatible schemas, so `vcat` /
`append!` work out of the box.

## Quick reference

| Task                              | Page                          |
| --------------------------------- | ----------------------------- |
| Plummer cluster                   | [Plummer](@ref)               |
| Galactic disk                     | [Exponential Disk](@ref)      |
| Galactic bulge                    | [Bulge](@ref)                 |
| SPH gas cloud                     | [Gas Cloud](@ref)             |
| Generic spherical system          | [Spherical System](@ref)      |
| Solar-system snapshot             | [Solar System](@ref)          |
| Translate / center / frame        | [Tools](@ref)                 |
| Observational data loaders        | [Data](@ref)                  |