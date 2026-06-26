# Types

All initial-condition generators are subtypes of `InitialConditionConfig`. Use
the [`generate`](@ref) dispatcher to produce a `StructArray` of particles
from any of these.

## Index

```@index
Pages = ["Types.md"]
```

## Initial-condition configurations

```@docs
PlummerStarCluster
ExponentialDisc
Bulge
GasCloud
SphericalSystem
```