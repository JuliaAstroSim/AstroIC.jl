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
        TotalMass::Number = 8.57u"Msun",
        ScaleRadius::Number = 0.075u"kpc",
        CutRadius::Number = 2.1u"kpc",
        q::Number = 0.5,
        α::Number = 1.8,
    )
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

"""
$(TYPEDSIGNATURES)

"""
function generate(config::Bulge, units = uAstro;
    RotationCurve = nothing,
    MaxRadius = 5 * config.ScaleRadius,
    MaxHeight = MaxRadius,
)
    uLen = getuLength(units)
    uVel = getuVel(units)

    NumSamples = config.NumSamples

    if MaxRadius < config.ScaleRadius
        @warn "`MaxRadius` is smaller than `ScaleRadius`, this may cause unphyiscal errors!"
    end

    R = eltype(config.ScaleRadius)[]
    z = eltype(config.ScaleRadius)[]
    
    while length(R) < NumSamples
        R_rand = randn() * config.CutRadius / sqrt(2)
        z_rand = randn() * config.q * config.CutRadius / sqrt(2)
        r_prime = sqrt(R_rand^2 + (z_rand/config.q)^2)
        if rand() <= (1 + r_prime / config.ScaleRadius)^(-config.α)
            push!(R, abs(R_rand))
            push!(z, z_rand)
        end
    end

    pos2d = StructArray(rand_pos_2d.(R))
    pos = StructArray(PVector.(pos2d.x, pos2d.y, z))

    #TODO: vel
    vel = [PVector(uVel) for i in 1:NumSamples]
    
    # Packing
    Promote = ismeasurement(config.ScaleRadius) || ismeasurement(config.TotalMass)
    Promote && @info "Promoting to `Measurement`"
    particles = StructArray(Star(units, id = i, Measurement=Promote) for i in 1:NumSamples)
    assign_particles(particles, :Pos, pos)
    assign_particles(particles, :Vel, vel)

    Mmean = config.TotalMass / NumSamples
    assign_particles(particles, :Mass, Mmean)

    return particles
end