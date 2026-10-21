"""
$(TYPEDEF)
$(TYPEDFIELDS)

`mass_shell` supports two forms:
- `Dict` or `DataFrame`: discrete radii `r` and shell mass per unit radius
  `m = dM/dr`. `Dict` keys may be either symbols (`:r`, `:m`) or strings
  (`"r"`, `"m"`).
- `Function`: an analytic radial function `r -> dM/dr`.

In both cases `m` must have units of mass per unit radius (e.g.
`Msun/kpc`), not a cumulative mass.
"""
struct SphericalSystem{I} <: InitialConditionConfig
    collection::Collection
    NumSamples::I

    mass_shell
end

function Base.show(io::IO, config::SphericalSystem)
    print(io,
        "Config of spherical system Initial Conditions:",
        "\n    Particle Collection: ", config.collection,
        "\n      Number of Samples: ", config.NumSamples,
    )
end

"""
$(TYPEDSIGNATURES)

`RotationCurve` is an optional `(radii, velocities)` tuple used to assign
rotational velocities. `VelocityDispersion` (a constant or a function
`r -> σ`) adds isotropic Gaussian velocity dispersion. Without either keyword
the generated velocities are zero.
"""
function generate(config::SphericalSystem, units = uAstro;
    MaxRadius = 0.0u"kpc",
    RotationCurve = nothing,
    VelocityDispersion = nothing,
    k = 2,
    bc = "nearest",
    rotational_ratio = 0.9,
    rng::AbstractRNG = Random.default_rng(),
)
    uLen = getuLength(units)
    uVel = getuVel(units)
    NumSamples = config.NumSamples

    if MaxRadius <= zero(MaxRadius)
        error("Keyword `MaxRadius` must be positive")
    end

    if config.mass_shell isa Function
        R = rejection_sampling(config.mass_shell, MaxRadius, NumSamples; rng = rng)
        x = collect(LinRange(MaxRadius/10000, MaxRadius, 10000))
        TotalMass = PhysicalParticles.NumericalIntegration.integrate(
            x, config.mass_shell.(x)
        )
    else
        # Discrete table: `m` is the shell mass per unit radius (dM/dr).
        if config.mass_shell isa Dict
            rkey = haskey(config.mass_shell, :r) ? :r : "r"
            mkey = haskey(config.mass_shell, :m) ? :m : "m"
            if !haskey(config.mass_shell, rkey) || !haskey(config.mass_shell, mkey)
                error("`mass_shell` of `SphericalSystem` must be `Dict` or `DataFrame` (with fields `:r`, `:m`) or `Function`")
            end
            r_table = config.mass_shell[rkey]
            m_table = config.mass_shell[mkey]
        elseif config.mass_shell isa DataFrame
            if !hasproperty(config.mass_shell, :r) || !hasproperty(config.mass_shell, :m)
                error("`mass_shell` of `SphericalSystem` must be `Dict` or `DataFrame` (with fields `:r`, `:m`) or `Function`")
            end
            r_table = config.mass_shell.r
            m_table = config.mass_shell.m
        else
            error("`mass_shell` of `SphericalSystem` must be `Dict` or `DataFrame` (with fields `:r`, `:m`) or `Function`")
        end

        # Evaluate the spline and sample in the table's length unit so that
        # the table, rejection box, and integrated mass share one unit.
        uR = unit(r_table[1])
        uM = unit(m_table[1])
        MaxRadius_table = uconvert(uR, MaxRadius)
        R = rejection_sampling(r_table, m_table, MaxRadius_table, NumSamples; rng = rng)

        x = collect(LinRange(MaxRadius_table/10000, MaxRadius_table, 10000))
        spl = Spline1D(ustrip.(uR, r_table), ustrip.(uM, m_table))
        TotalMass = PhysicalParticles.NumericalIntegration.integrate(
            x, spl(ustrip.(uR, x)) * uM
        )
    end

    pos = [rand_pos_3d(ri; rng = rng) for ri in R]

    vel = [PVector(uVel) for i in 1:NumSamples]
    if !isnothing(RotationCurve)
        xc, vc = RotationCurve
        spl = Spline1D(ustrip.(uLen, xc), ustrip.(uVel, vc); k, bc)
        v = spl(ustrip.(uLen, R)) * uVel
        px = [p.x for p in pos]
        py = [p.y for p in pos]
        vel = [rotational_velocity(x, y, vv, rotational_ratio; rng = rng)
               for (x, y, vv) in zip(px, py, v)]
    end

    if !isnothing(VelocityDispersion)
        for i in eachindex(vel)
            σ = _velocity_dispersion(VelocityDispersion, norm(pos[i]))
            vel[i] = vel[i] + PVector(randn(rng) * σ, randn(rng) * σ, randn(rng) * σ)
        end
    end

    if !isnothing(RotationCurve) || !isnothing(VelocityDispersion)
        # Cancel out any net drift
        v0 = mean(vel)
        vel = vel .- v0
        vel = uconvert.(uVel, vel)
    end

    # Packing
    particles = StructArray(Star(units, id = i, collection = config.collection) for i in 1:NumSamples)
    assign_particles(particles, :Pos, pos)
    assign_particles(particles, :Vel, vel)

    Mmean = TotalMass / NumSamples
    assign_particles(particles, :Mass, Mmean)

    return particles
end