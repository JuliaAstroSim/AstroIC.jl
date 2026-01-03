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

"""
    generate_milkyway_baryon_particles(Np)

Generate Milky Way baryon particles (bulge, thin/thick discs, HI/HII gas).
"""
function generate_milkyway_baryon_particles(Np)
    TotalMass_bulge = 8.5708e9u"Msun"
    TotalMass_thin = uconvert(u"Msun", 2pi * 1003.12u"Msun/pc^2" * (2.42u"kpc")^2) # 3.691165259383738e10 M⊙
    TotalMass_thick = uconvert(u"Msun", 2pi * 167.93u"Msun/pc^2" * (3.17u"kpc")^2) # 1.0602949202938915e10 M⊙
    TotalMass_HI = 1.0674e10u"Msun"
    TotalMass_HII = 1.2303e9u"Msun"
    TotalMass_baryons = TotalMass_bulge + TotalMass_thin + TotalMass_thick + TotalMass_HI + TotalMass_HII

    NumSamples_bulge = ceil(Int, TotalMass_bulge/TotalMass_baryons * Np)
    NumSamples_thin  = ceil(Int, TotalMass_thin/TotalMass_baryons * Np)
    NumSamples_thick = ceil(Int, TotalMass_thick/TotalMass_baryons * Np)
    NumSamples_HI    = ceil(Int, TotalMass_HI/TotalMass_baryons * Np)
    NumSamples_HII = Np - NumSamples_bulge - NumSamples_thin - NumSamples_thick - NumSamples_HI

    if iszero(NumSamples_HII)
        error("Number of HII particles is zero. Try increase number of particles")
    end

    @info "NumSamples of bulge: $(NumSamples_bulge)"
    @info "NumSamples of thin:  $(NumSamples_thin)"
    @info "NumSamples of thick: $(NumSamples_thick)"
    @info "NumSamples of HI:    $(NumSamples_HI)"
    @info "NumSamples of HII:   $(NumSamples_HII)"

    particles_bulge = generate(Bulge(;
        collection = STAR,
        NumSamples = NumSamples_bulge,
        TotalMass = TotalMass_bulge,
        ScaleRadius = 0.075u"kpc",
        CutRadius = 2.1u"kpc",
        q = 0.5,
        α = 1.8,
    ))

    particles_stellar_thin = generate(ExponentialDisc(;
        collection = STAR,
        NumSamples = NumSamples_thin,
        TotalMass = TotalMass_thin,
        ScaleRadius = 2.42u"kpc",
        ScaleHeight = 0.3u"kpc",
    ); RotationCurve = nothing)

    particles_stellar_thick = generate(ExponentialDisc(;
        collection = STAR,
        NumSamples = NumSamples_thick,
        TotalMass = TotalMass_thick,
        ScaleRadius = 3.17u"kpc",
        ScaleHeight = 0.9u"kpc",
    ); RotationCurve = nothing)

    particles_gas_HI = generate(ExponentialDisc(;
        collection = STAR,
        NumSamples = NumSamples_HI,
        TotalMass = TotalMass_HI,
        ScaleRadius = 7.0u"kpc",
        ScaleHeight = 0.085u"kpc",
        HoleRadius = 4.0u"kpc",
    ))

    particles_gas_HII = generate(ExponentialDisc(;
        collection = STAR,
        NumSamples = NumSamples_HII,
        TotalMass = TotalMass_HII,
        ScaleRadius = 1.5u"kpc",
        ScaleHeight = 0.045u"kpc",
        HoleRadius = 12.0u"kpc",
    ))

    particles = [particles_bulge; particles_stellar_thin; particles_stellar_thick; particles_gas_HI; particles_gas_HII]
    return particles
end