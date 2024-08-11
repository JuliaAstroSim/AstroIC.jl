"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
mutable struct ExponentialDisk{I, Len, MASS} <: InitialConditionConfig
    collection::Collection
    NumSamples::I

    TotalMass::MASS

    ScaleRadius::Len
    ScaleHeight::Len
    HoleRadius::Len
end

"""
$(TYPEDSIGNATURES)

- `collection` particle type
- `NumSamples` amount of particles
"""
function ExponentialDisk(;
        collection::Collection = STAR,
        NumSamples::Int64 = 1000,
        TotalMass::Number = 1.0e10u"Msun",
        ScaleRadius::Number = 2.0u"kpc",
        ScaleHeight::Number = 0.02u"kpc",
        HoleRadius::Number = 0.0u"kpc",
    )
    return ExponentialDisk(collection, NumSamples, TotalMass, ScaleRadius, ScaleHeight, HoleRadius)
end

function Base.show(io::IO, config::ExponentialDisk)
    print(io,
        "Config of Exponential Disk Initial Conditions:",
        "\n    Particle Collection: ", config.collection,
        "\n      Number of Samples: ", config.NumSamples,
        "\n           Scale Radius: ", config.ScaleRadius,
        "\n           Scale Height: ", config.ScaleHeight,
        "\n            Hole Radius: ", config.HoleRadius,
        "\n             Total Mass: ", config.TotalMass,
    )
end

"""
$(TYPEDSIGNATURES)

Ref: https://en.wikipedia.org/wiki/Exponential_distribution
"""
function exponential_pdf(beta, x)
    if x > zero(x)
        y = exp(-x/beta)/beta
    else
        y = zero(x)
    end
    return y
end

"""
$(TYPEDSIGNATURES)

Ref: https://en.wikipedia.org/wiki/Exponential_distribution
"""
function exponential_cdf(beta, x)
    if x > zero(x)
        y = 1-exp(-x/beta)
    else
        y = zero(x)
    end
    return y
end

"""
$(TYPEDSIGNATURES)

Ref: https://en.wikipedia.org/wiki/Exponential_distribution
"""
function exponential_cdf_inv(beta, y)
    x = -beta*log(1-y)
end

"""
$(TYPEDSIGNATURES)

"""
function exponential_pdf_hole(beta, Rm, x)
    if x > zero(x)
        y = exp(-x/beta-Rm/x)/beta
    else
        y = zero(x)
    end
    return y
end

"""
$(TYPEDSIGNATURES)

"""
function sech2_pdf(x)
    return sech(x)^2
end

"""
$(TYPEDSIGNATURES)

"""
function sech2_cdf(x)
    return 0.5*tanh(x)+0.5
end

"""
$(TYPEDSIGNATURES)

"""
function sech2_cdf_inv(u)
    return atanh(2*u-1)
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

    # generate radii and heights
    if MaxRadius < config.ScaleRadius
        @warn "`MaxRadius` is smaller than `ScaleRadius`, this may cause unphyiscal errors!"
    end
    
    if MaxHeight < config.ScaleHeight
        @warn "`MaxHeight` is smaller than `ScaleHeight`, this may cause unphyiscal errors!"
    end

    if iszero(config.HoleRadius)
        R = exponential_cdf_inv.(config.ScaleRadius, rand(NumSamples))
        for i in eachindex(R) # filter with MaxRadius, make sure that radii are all inbound
            while R[i] > MaxRadius
                R[i] = exponential_cdf_inv(config.ScaleRadius, rand())
            end
        end

        z = exponential_cdf_inv.(config.ScaleHeight, rand(NumSamples)) .* rand([-1,1], NumSamples)
        for i in eachindex(z)
            while z[i] > MaxHeight
                z[i] = exponential_cdf_inv(config.ScaleHeight, rand())
            end
        end
    else
        @info "Generating exponential disk with hole"
        R = eltype(config.ScaleRadius)[]
        while length(R) < NumSamples
            x = exponential_cdf_inv(config.ScaleRadius, rand())
            if rand() <= exponential_pdf_hole(config.ScaleRadius, config.HoleRadius, x) / exponential_pdf(config.ScaleRadius, x) # Reject sampling
                push!(R, x)
            end
        end

        # sech² Hyperbolic secant-squared distribution
        z = sech2_cdf_inv.(rand(NumSamples)) * 2 * config.ScaleHeight
        for i in eachindex(z)
            while z[i] > MaxHeight
                z[i] = sech2_cdf_inv(rand()) * 2 * config.ScaleHeight
            end
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