"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
mutable struct Bulge{I, Len, MASS} <: InitialConditionConfig
    collection::Collection
    NumSamples::I

    TotalMass::MASS

    ScaleRadius::Len
    CutRadius::Len
    q
    α
end

"""
$(TYPEDSIGNATURES)

- `collection` particle type
- `NumSamples` amount of particles
"""
function Bulge(;
        collection::Collection = STAR,
        NumSamples::Int64 = 1000,
        TotalMass::Number = 8.57e9u"Msun",
        ScaleRadius::Number = 0.075u"kpc",
        CutRadius::Number = 2.1u"kpc",
        q::Number = 0.5,
        α::Number = 1.8,
    )
    # `ScaleRadius` and `CutRadius` share one type parameter; normalize them
    # so callers may pass either radius in any length unit.
    CutRadius = uconvert(unit(ScaleRadius), CutRadius)
    return Bulge(collection, NumSamples, TotalMass, ScaleRadius, CutRadius, q, α)
end

function Base.show(io::IO, config::Bulge)
    print(io,
        "Config of Bulge Initial Conditions:",
        "\n    Particle Collection: ", config.collection,
        "\n      Number of Samples: ", config.NumSamples,
        "\n           Scale Radius: ", config.ScaleRadius,
        "\n             Cut Radius: ", config.CutRadius,
        "\n             Total Mass: ", config.TotalMass,
        "\n         Axis Radio (q): ", config.q,
        "\n                      α: ", config.α,
    )
end

function pdf(config::Bulge, R, z)
    r_prime = sqrt(R^2 + (z / config.q)^2)
    return exp(-r_prime^2/ustrip(config.CutRadius)^2) / (1 + r_prime/ustrip(config.ScaleRadius))^config.α * 2π*R
end

"""
$(TYPEDSIGNATURES)

`MaxRadius` and `MaxHeight` truncate the density profile. `TotalMass` is
distributed over all sampled particles, so the truncated density is
renormalized to the requested total mass.

Passing `RotationCurve = (radii, velocities)` assigns rotational velocities.
`VelocityDispersion` (a constant or a function `r -> σ`) adds isotropic
Gaussian velocity dispersion. Without either keyword the velocities are zero.
"""
function generate(config::Bulge, units = uAstro;
    RotationCurve = nothing,
    VelocityDispersion = nothing,
    MaxRadius = 5 * config.ScaleRadius,
    MaxHeight = MaxRadius,
    k = 2,
    bc = "nearest",
    rotational_ratio = 0.9,
    rng::AbstractRNG = Random.default_rng(),
)
    uLen = getuLength(units)
    uVel = getuVel(units)

    # Sample in the config's length unit; convert the bounds to it so the
    # dimensionless rejection-sampling numbers and `pdf` agree.
    uScale = unit(config.ScaleRadius)
    MaxRadius = uconvert(uScale, MaxRadius)
    MaxHeight = uconvert(uScale, MaxHeight)

    NumSamples = config.NumSamples

    if MaxRadius <= zero(MaxRadius)
        error("Keyword `MaxRadius` must be positive")
    end
    if MaxHeight <= zero(MaxHeight)
        error("Keyword `MaxHeight` must be positive")
    end

    if MaxRadius < config.ScaleRadius
        @warn "`MaxRadius` is smaller than `ScaleRadius`, this may cause unphyiscal errors!"
    end

    R = eltype(config.ScaleRadius)[]
    sizehint!(R, NumSamples)
    z = eltype(config.ScaleRadius)[]
    sizehint!(z, NumSamples)

    target(xy) = -pdf(config,xy[1],xy[2])
    pdf_maximum = -minimum_func(target, [ustrip(config.ScaleRadius), ustrip(config.ScaleRadius)*config.q])[1]

    # rejection sampling
    while length(R) < NumSamples
        R_rand = rand(rng) * ustrip(uScale, MaxRadius)
        z_rand = rand(rng) * 2*ustrip(uScale, MaxHeight) - ustrip(uScale, MaxHeight)
        if rand(rng) < pdf(config, R_rand, z_rand) / pdf_maximum
            push!(R, R_rand * uScale)
            push!(z, z_rand * uScale)
        end
    end

    # while length(R) < NumSamples
    #     R_rand = randn() * config.CutRadius / sqrt(2)
    #     z_rand = randn() * config.q * config.CutRadius / sqrt(2)
    #     r_prime = sqrt(R_rand^2 + (z_rand/config.q)^2)
    #     if rand() <= (1 + r_prime / config.ScaleRadius)^(-config.α)
    #         push!(R, abs(R_rand))
    #         push!(z, z_rand)
    #     end
    # end

    pos2d = StructArray([rand_pos_2d(ri; rng = rng) for ri in R])
    pos = StructArray(PVector.(pos2d.x, pos2d.y, z))

    vel = [PVector(uVel) for i in 1:NumSamples]
    if !isnothing(RotationCurve)
        xc, vc = RotationCurve
        spl = Spline1D(ustrip.(uLen, xc), ustrip.(uVel, vc); k, bc)
        v = spl(ustrip.(uLen, R)) * uVel
        vel = [rotational_velocity(px, py, vv, rotational_ratio; rng = rng)
               for (px, py, vv) in zip(pos.x, pos.y, v)]
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
    Promote = ismeasurement(config.ScaleRadius) || ismeasurement(config.TotalMass)
    Promote && @info "Promoting to `Measurement`"
    particles = StructArray(Star(units, id = i, Measurement=Promote, collection = config.collection) for i in 1:NumSamples)
    assign_particles(particles, :Pos, pos)
    assign_particles(particles, :Vel, vel)

    Mmean = config.TotalMass / NumSamples
    assign_particles(particles, :Mass, Mmean)

    return particles
end