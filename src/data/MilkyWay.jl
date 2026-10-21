"""
$(TYPEDSIGNATURES)
"""
function load_MW_RC_Eilers2019(filename = joinpath(@__DIR__, "MilkyWay/MW_RC_Eilers2019.dat"))
    dfEilers2019 = DataFrame(CSV.File(filename; header = false, delim=" ", ignorerepeated = true))
    rename!(dfEilers2019, [:r, :v, :σ_low, :σ_high])
    return dfEilers2019
end

"""
$(TYPEDSIGNATURES)
"""
function load_MW_RC_Mroz2019(filename = joinpath(@__DIR__, "MilkyWay/MW_RC_Mroz2019.dat"))
    dfMroz2019 = DataFrame(CSV.File(filename; header = false, delim=" ", ignorerepeated = true))
    rename!(dfMroz2019, [:r, :v, :σ_low, :σ_high])
    return dfMroz2019
end

"""
$(TYPEDSIGNATURES)
"""
function load_MW_RC_stddev_W21(filename = joinpath(@__DIR__, "MilkyWay/MW_RC_stddev_W21.txt"))
    dfstddev_W21 = DataFrame(CSV.File(filename; header = false, delim=" ", ignorerepeated = true))
    rename!(dfstddev_W21, [:id, :r, :v, :σ_v, :σ_r])
    return dfstddev_W21
end

"""
$(TYPEDSIGNATURES)
"""
function load_MW_RC_DS_W21(filename = joinpath(@__DIR__, "MilkyWay/MW_dsvel-W21-RC.txt"))
    dfDS_W21 = DataFrame(CSV.File(filename; header = false, skipto=3, delim="\t", ignorerepeated = true))
    rename!(dfDS_W21, [:r, :v_b, :v_CDM, :v_QUMOND, :v_MOG])
    return dfDS_W21
end

function _mw_default_rotation_curve()
    df = load_MW_RC_Eilers2019()
    return (df.r .* u"kpc", df.v .* u"km/s")
end

"""
    generate_milkyway_baryon_particles(Np; RotationCurve = _mw_default_rotation_curve())

Generate Milky Way baryon particles (bulge, thin/thick discs, HI/H2 gas).
The mass model and gas-disc parameters follow McMillan (2017, MNRAS 465, 76).

`RotationCurve` is a `(radii, velocities)` tuple applied to the stellar discs
and gas components. By default the Eilers et al. (2019) rotation curve bundled
with the package is used; pass `nothing` to generate non-rotating discs/gas.
The bulge is dispersion-supported with a constant 100 km/s (1D) velocity
dispersion. `rng` is forwarded to every component generator.
"""
function generate_milkyway_baryon_particles(Np;
        RotationCurve = _mw_default_rotation_curve(),
        rng::AbstractRNG = Random.default_rng())
    TotalMass_bulge = 8.5708e9u"Msun"
    TotalMass_thin = uconvert(u"Msun", 2pi * 1003.12u"Msun/pc^2" * (2.42u"kpc")^2) # 3.691165259383738e10 M⊙
    TotalMass_thick = uconvert(u"Msun", 2pi * 167.93u"Msun/pc^2" * (3.17u"kpc")^2) # 1.0602949202938915e10 M⊙
    TotalMass_HI = 1.0674e10u"Msun"
    TotalMass_H2 = 1.2303e9u"Msun"
    TotalMass_baryons = TotalMass_bulge + TotalMass_thin + TotalMass_thick + TotalMass_HI + TotalMass_H2

    NumSamples_bulge = ceil(Int, TotalMass_bulge/TotalMass_baryons * Np)
    NumSamples_thin  = ceil(Int, TotalMass_thin/TotalMass_baryons * Np)
    NumSamples_thick = ceil(Int, TotalMass_thick/TotalMass_baryons * Np)
    NumSamples_HI    = ceil(Int, TotalMass_HI/TotalMass_baryons * Np)
    NumSamples_H2    = Np - NumSamples_bulge - NumSamples_thin - NumSamples_thick - NumSamples_HI

    if NumSamples_H2 < 1
        error("Number of H2 particles is $(NumSamples_H2). Increase the total number of particles")
    end

    @info "NumSamples of bulge: $(NumSamples_bulge)"
    @info "NumSamples of thin:  $(NumSamples_thin)"
    @info "NumSamples of thick: $(NumSamples_thick)"
    @info "NumSamples of HI:    $(NumSamples_HI)"
    @info "NumSamples of H2:    $(NumSamples_H2)"

    particles_bulge = generate(Bulge(;
        collection = STAR,
        NumSamples = NumSamples_bulge,
        TotalMass = TotalMass_bulge,
        ScaleRadius = 0.075u"kpc",
        CutRadius = 2.1u"kpc",
        q = 0.5,
        α = 1.8,
    ); MaxRadius = 5*2.1u"kpc", VelocityDispersion = 100.0u"km/s", rng = rng)

    particles_stellar_thin = generate(ExponentialDisc(;
        collection = STAR,
        NumSamples = NumSamples_thin,
        TotalMass = TotalMass_thin,
        ScaleRadius = 2.42u"kpc",
        ScaleHeight = 0.3u"kpc",
    ); RotationCurve = RotationCurve, rng = rng)

    particles_stellar_thick = generate(ExponentialDisc(;
        collection = STAR,
        NumSamples = NumSamples_thick,
        TotalMass = TotalMass_thick,
        ScaleRadius = 3.17u"kpc",
        ScaleHeight = 0.9u"kpc",
    ); RotationCurve = RotationCurve, rng = rng)

    particles_gas_HI = generate(ExponentialDisc(;
        collection = GAS,
        NumSamples = NumSamples_HI,
        TotalMass = TotalMass_HI,
        ScaleRadius = 7.0u"kpc",
        ScaleHeight = 0.085u"kpc",
        HoleRadius = 4.0u"kpc",
    ); RotationCurve = RotationCurve, rng = rng)

    particles_gas_H2 = generate(ExponentialDisc(;
        collection = GAS,
        NumSamples = NumSamples_H2,
        TotalMass = TotalMass_H2,
        ScaleRadius = 1.5u"kpc",
        ScaleHeight = 0.045u"kpc",
        HoleRadius = 12.0u"kpc",
    ); RotationCurve = RotationCurve, rng = rng)

    particles = [particles_bulge; particles_stellar_thin; particles_stellar_thick; particles_gas_HI; particles_gas_H2]
    return particles
end