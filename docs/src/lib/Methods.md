# Methods

## Index

```@index
Pages = ["Methods.md"]
```

## Initial-condition generation

```@docs
generate
generate(::PlummerStarCluster)
generate(::ExponentialDisc)
generate(::Bulge)
generate(::GasCloud)
generate(::SphericalSystem)
solarsystem
```

## Particle-system manipulation

```@docs
setpos
setvel
addpos
addvel
```

## Geometry & physics helpers

```@docs
gridpoints
helio2xyz
vmean
vmean2
```

## Data loaders

```@docs
load_data_MW_satellites
load_data_UFDs
load_SPARC_LTGs_RC
load_SPARC_LTGs_data
load_li2018_SPARC
load_SPARC_Xray_ETGs_data
load_SPARC_rotating_ETGs_data
load_SPARC_rotating_ETGs_RC
load_SPARC_rotating_ETGs_rotmod
load_MW_RC_Eilers2019
load_MW_RC_Mroz2019
load_MW_RC_stddev_W21
load_MW_RC_DS_W21
load_massive_dwarf_CO_RC
load_massive_dwarf_DM_RC
load_massive_dwarf_Baryon_RC
generate_milkyway_baryon_particles
```

## Internal helpers

These functions are not part of the public API but are used internally by the
generators above. They are listed here so that Documenter does not complain
about undocumented docstrings, and so that interested users can find them.

```@autodocs
Modules = [AstroIC]
Public = false
Filter = t -> nameof(t) in (
    :rand_pos_2d,
    :rand_pos_3d,
    :rotational_velocity,
    :rotational_velocity_acc,
    :freefall_velocity_acc,
    :minimum_func,
    :exponential_pdf,
    :exponential_cdf,
    :exponential_cdf_inv,
    :exponential_pdf_hole,
    :sech2_pdf,
    :sech2_cdf,
    :sech2_cdf_inv,
    :rejection_sampling,
    :Rd_from_RHI,
)
```