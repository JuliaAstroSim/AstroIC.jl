# Data

AstroIC ships with CSV / DAT files from a number of observational surveys so
that you can drive simulations directly from published measurements.

## Milky Way satellite galaxies

```@repl data_mw
using AstroIC
df = load_data_MW_satellites()
size(df)
```

See [`load_data_MW_satellites`](@ref), [`load_data_UFDs`](@ref) and
[`generate_milkyway_baryon_particles`](@ref) in the
[Library reference](@ref Methods) for full signatures.

`generate_milkyway_baryon_particles(Np)` produces a composite Milky-Way-like
system of bulge + thin disc + thick disc + HI + HII gas using the published
mass fractions, sampled from the corresponding `Bulge` / `ExponentialDisc`
configs.

## SPARC rotation curves

The Spitzer Photometry and Accurate Rotation Curves (SPARC, [Lelli+ 2016](https://arxiv.org/abs/1606.09251))
database covers 175 late-type galaxies.

```@repl data_sparc
using AstroIC

df_sparc = load_SPARC_LTGs_data()        # full galaxy parameter table
df_rc    = load_SPARC_LTGs_RC()          # mass-model rotation curves
df_li18  = load_li2018_SPARC()           # Radial Acceleration Relation
```

For early-type galaxies (Lelli+ 2017):

```@repl data_sparc_etg
using AstroIC

df_xray = load_SPARC_Xray_ETGs_data()
df_etg  = load_SPARC_rotating_ETGs_data()

# Per-galaxy rotation curves (bulge / disk decomposition)
rc_bulge = load_SPARC_rotating_ETGs_RC("NGC2685", :bulge)
rc_disk  = load_SPARC_rotating_ETGs_RC("NGC2824", :disk)
rotmod   = load_SPARC_rotating_ETGs_rotmod("NGC2685")
```

See [`load_SPARC_LTGs_RC`](@ref), [`load_SPARC_LTGs_data`](@ref),
[`load_li2018_SPARC`](@ref), [`load_SPARC_Xray_ETGs_data`](@ref),
[`load_SPARC_rotating_ETGs_data`](@ref),
[`load_SPARC_rotating_ETGs_RC`](@ref) and
[`load_SPARC_rotating_ETGs_rotmod`](@ref) in the
[Library reference](@ref Methods) for full signatures.

## Milky Way rotation curves

Three independent measurements of the Milky Way rotation curve are bundled:

| Function                  | Reference                              |
| ------------------------- | -------------------------------------- |
| `load_MW_RC_Eilers2019`   | Eilers+ (2019), [ApJ 871, 120](https://iopscience.iop.org/article/10.3847/1538-4357/aaf648)         |
| `load_MW_RC_Mroz2019`     | Mróz+ (2019), [ApJ 870, L10](https://iopscience.iop.org/article/10.3847/2041-8213/aaf73f)           |
| `load_MW_RC_DS_W21`       | Wang+ (2021), comparing CDM / QUMOND / MOG fits |
| `load_MW_RC_stddev_W21`   | Wang+ (2021) standard-deviation table           |

```julia
using AstroIC
df_eilers = load_MW_RC_Eilers2019()
df_mroz   = load_MW_RC_Mroz2019()
df_w21    = load_MW_RC_DS_W21()
df_w21σ   = load_MW_RC_stddev_W21()
```

See [`load_MW_RC_Eilers2019`](@ref), [`load_MW_RC_Mroz2019`](@ref),
[`load_MW_RC_stddev_W21`](@ref) and [`load_MW_RC_DS_W21`](@ref) in the
[Library reference](@ref Methods) for full signatures.

## Massive dwarf galaxy rotation curves (Cooke 2022)

Six massive dwarf rotation-curve decompositions are bundled
([Cooke 2022, MNRAS 512, 1012](https://arxiv.org/abs/2203.00694)):

```julia
using AstroIC
co  = load_massive_dwarf_CO_RC("NGC1035")        # CO gas contribution
dm  = load_massive_dwarf_DM_RC("NGC1035")        # dark matter
baryon = load_massive_dwarf_Baryon_RC("NGC1035") # total baryonic
```

See [`load_massive_dwarf_CO_RC`](@ref), [`load_massive_dwarf_DM_RC`](@ref) and
[`load_massive_dwarf_Baryon_RC`](@ref) in the
[Library reference](@ref Methods) for full signatures.