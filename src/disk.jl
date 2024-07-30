mutable struct ExponentialDisk{I, Len, MASS} <: InitialConditionConfig
    collection::Collection
    NumSamples::I

    TotalMass::MASS

    ScaleRadius::Len
    ScaleHeight::Len
end

"""
    struct ExponentialDisk

## Fields

- `collection` particle type
- `NumSamples` amount of particles
"""
function ExponentialDisk(;
        collection::Collection = STAR,
        NumSamples::Int64 = 1000,
        TotalMass::Number = 1.0e10u"Msun",
        ScaleRadius::Number = 2.0u"kpc",
        ScaleHeight::Number = 0.02u"kpc",
    )
    return ExponentialDisk(collection, NumSamples, TotalMass, ScaleRadius, ScaleHeight)
end

function Base.show(io::IO, config::ExponentialDisk)
    print(io,
        "Config of Exponential Disk Initial Conditions:",
        "\n    Particle Collection: ", config.collection,
        "\n      Number of Samples: ", config.NumSamples,
        "\n           Scale Radius: ", config.ScaleRadius,
        "\n           Scale Height: ", config.ScaleHeight,
        "\n             Total Mass: ", config.TotalMass,
    )
end


"""
$(TYPEDSIGNATURES)

Ref: https://en.wikipedia.org/wiki/Exponential_distribution
"""
function exponential_pdf_inv(beta, y)
    x = -beta*log(1-y)
end

"""
$(TYPEDSIGNATURES)

"""
function generate(config::ExponentialDisk, units = uAstro;
        RotationCurve = nothing,
        MaxRadius = 5 * config.ScaleRadius,
        MaxHeight = MaxRadius,
        k = 1, # degree of rotation curve interpolation/extrapolation spline (1 = linear, 2 = quadratic, 3 = cubic, up to 5)
        bc = "nearest", # behavior when evaluating the spline outside the support domain, which is (minimum(x), maximum(x)). The allowed values are "nearest", "zero", "extrapolate", "error"
        rotational_ratio = 0,
    )
    uLen = getuLength(units)
    uVel = getuVel(units)

    NumSamples = config.NumSamples
    ScaleRadius = uconvert(uLen, config.ScaleRadius)
    ScaleHeight = uconvert(uLen, config.ScaleHeight)

    # generate radii and heights
    if MaxRadius < ScaleRadius
        @warn "`MaxRadius` is smaller than `ScaleRadius`, this may cause unphyiscal errors!"
    end

    R = exponential_pdf_inv.(ScaleRadius, rand(NumSamples))
    for i in eachindex(R)
        while R[i] > MaxRadius
            R[i] = exponential_pdf_inv(ScaleRadius, rand())
        end
    end

    if MaxHeight < ScaleHeight
        @warn "`MaxHeight` is smaller than `ScaleHeight`, this may cause unphyiscal errors!"
    end

    z = exponential_pdf_inv.(ScaleHeight, rand(NumSamples))
    for i in eachindex(z)
        while z[i] > MaxHeight
            z[i] = exponential_pdf_inv(ScaleHeight, rand())
        end
    end

    pos2d = StructArray(rand_pos_2d.(R))
    pos = StructArray(PVector.(pos2d.x, pos2d.y, z))

    # Extrapolate rotation curve
    if isnothing(RotationCurve)
        vel = [PVector(uVel) for i in 1:NumSamples]
    else
        xc, vc = RotationCurve
        spl = Spline1D(ustrip.(unit(eltype(xc)), xc), vc; k, bc)
        v = spl(ustrip.(unit(eltype(R)), R))
        vel = rotational_velocity.(pos.x, pos.y, pos.z, v, rotational_ratio)
    end

    # Packing
    Promote = ismeasurement(config.ScaleRadius) || ismeasurement(config.ScaleHeight) || ismeasurement(config.TotalMass)
    Promote && @info "Promoting to `Measurement`"
    particles = StructArray(Star(units, id = i, Measurement=Promote) for i in 1:NumSamples)
    assign_particles(particles, :Pos, pos)
    assign_particles(particles, :Vel, vel)

    Mmean = config.TotalMass / NumSamples
    assign_particles(particles, :Mass, Mmean)

    return particles
end