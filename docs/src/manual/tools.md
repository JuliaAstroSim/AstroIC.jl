# Tools

The `Tools` module provides helpers for translating, centering, and re-framing
a particle system after generation.

See [`addpos`](@ref), [`addvel`](@ref), [`setpos`](@ref) and [`setvel`](@ref)
in the [Library reference](@ref Methods) for full signatures.

## Example: colliding Plummer spheres

A typical N-body initial-condition is several Plummer clusters placed at
different positions with different bulk velocities, concatenated into a single
particle array:

```@repl tools
using AstroIC
using PhysicalParticles, UnitfulAstro

config = PlummerStarCluster(NumSamples = 1000)
galaxy1 = generate(config);
galaxy2 = generate(config);
galaxy3 = generate(config);
galaxy4 = generate(config);

setpos(galaxy1, PVector(-0.1, 0.0, 0.0, u"kpc"))
setpos(galaxy2, PVector(+0.1, 0.0, 0.0, u"kpc"))
setpos(galaxy3, PVector(0.0, +0.1, 0.0, u"kpc"))
setpos(galaxy4, PVector(0.0, -0.1, 0.0, u"kpc"))

setvel(galaxy1, PVector(+0.4, +0.8, 0.0, u"kpc/Gyr"))
setvel(galaxy2, PVector(-0.4, -1.0, -0.2, u"kpc/Gyr"))
setvel(galaxy3, PVector(+1.0, -0.4, 0.4, u"kpc/Gyr"))
setvel(galaxy4, PVector(-1.0, +0.4, -0.2, u"kpc/Gyr"))

data = deepcopy(galaxy1);
append!(data, galaxy2);
append!(data, galaxy3);
append!(data, galaxy4);
data.ID .= 1:4000;

length(data)
```

## When to use which

- [`setpos`](@ref) / [`setvel`](@ref) — place the system so its **median**
  position (or mass-weighted mean velocity) is at a target value. Useful when
  the original sampling already scatters particles around zero.
- [`addpos`](@ref) / [`addvel`](@ref) — translate by a fixed vector,
  irrespective of where the system currently sits.

`setpos` uses `median(data, :Pos)` (robust against outliers), `setvel` uses
`averagebymass(data, :Vel)` (physical bulk velocity).