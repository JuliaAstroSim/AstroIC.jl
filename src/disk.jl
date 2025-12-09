"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
mutable struct ExponentialDisc{I, Len, MASS} <: InitialConditionConfig
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
function ExponentialDisc(;
        collection::Collection = STAR,
        NumSamples::Int64 = 1000,
        TotalMass::Number = 1.0e10u"Msun",
        ScaleRadius::Number = 2.0u"kpc",
        ScaleHeight::Number = 0.02u"kpc",
        HoleRadius::Number = 0.0u"kpc",
    )
    return ExponentialDisc(collection, NumSamples, TotalMass, ScaleRadius, ScaleHeight, HoleRadius)
end

function Base.show(io::IO, config::ExponentialDisc)
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

function pdf(config::ExponentialDisc, R, z)
    if iszero(config.HoleRadius)
        return exp(-abs(z)/ustrip(config.ScaleHeight) - R/ustrip(config.ScaleRadius)) / ustrip(config.ScaleHeight) * 2π * R
    else
        return exp(-ustrip(config.HoleRadius)/R - R/ustrip(config.ScaleRadius)) * sech(z/ustrip(config.ScaleHeight)/2)^2 * 2π * R
    end
end

"""
$(TYPEDSIGNATURES)

"""
function generate(config::ExponentialDisc, units = uAstro;
    RotationCurve = nothing,
    MaxRadius = 5 * config.ScaleRadius,
    MaxHeight = MaxRadius,
    k = 2, # degree of rotation curve interpolation/extrapolation spline (1 = linear, 2 = quadratic, 3 = cubic, up to 5)
    bc = "nearest", # behavior when evaluating the spline outside the support domain, which is (minimum(x), maximum(x)). The allowed values are "nearest", "zero", "extrapolate", "error"
    rotational_ratio = 0.9,
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

    # if iszero(config.HoleRadius)
    #     R = exponential_cdf_inv.(config.ScaleRadius, rand(NumSamples))
    #     for i in eachindex(R) # filter with MaxRadius, make sure that radii are all inbound
    #         while R[i] > MaxRadius
    #             R[i] = exponential_cdf_inv(config.ScaleRadius, rand())
    #         end
    #     end

    #     z = exponential_cdf_inv.(config.ScaleHeight, rand(NumSamples)) .* rand([-1,1], NumSamples)
    #     for i in eachindex(z)
    #         while z[i] > MaxHeight
    #             z[i] = exponential_cdf_inv(config.ScaleHeight, rand())
    #         end
    #     end
    # else
    #     # @info "Generating exponential disk with hole"
    #     R = eltype(config.ScaleRadius)[]
    #     while length(R) < NumSamples
    #         x = exponential_cdf_inv(config.ScaleRadius, rand())
    #         if rand() <= exponential_pdf_hole(config.ScaleRadius, config.HoleRadius, x) / exponential_pdf(config.ScaleRadius, x) # Reject sampling
    #             push!(R, x)
    #         end
    #     end

    #     # sech² Hyperbolic secant-squared distribution
    #     z = sech2_cdf_inv.(rand(NumSamples)) * 2 * config.ScaleHeight
    #     for i in eachindex(z)
    #         while z[i] > MaxHeight
    #             z[i] = sech2_cdf_inv(rand()) * 2 * config.ScaleHeight
    #         end
    #     end
    # end

    R = eltype(config.ScaleRadius)[]
    sizehint!(R, NumSamples)
    z = eltype(config.ScaleRadius)[]
    sizehint!(z, NumSamples)

    target(xy) = -pdf(config,xy[1],xy[2])
    if iszero(config.HoleRadius)
        pdf_maximum = -minimum_func(target, [ustrip(config.ScaleRadius), ustrip(config.ScaleRadius)])[1]
    else
        pdf_maximum = pdf(config, ustrip(sqrt(config.ScaleRadius * config.HoleRadius)), 0.0)
    end
    # @show pdf_maximum

    # rejection sampling
    while length(R) < NumSamples
        R_rand = rand() * ustrip(MaxRadius)
        z_rand = rand() * 2*ustrip(MaxHeight) - ustrip(MaxHeight)
        if rand() < pdf(config, R_rand, z_rand) / pdf_maximum
            push!(R, R_rand * uLen)
            push!(z, z_rand * uLen)
        end
    end

    pos2d = StructArray(rand_pos_2d.(R))
    pos = StructArray(PVector.(pos2d.x, pos2d.y, z))

    # Extrapolate rotation curve
    if isnothing(RotationCurve)
        vel = [PVector(uVel) for i in 1:NumSamples]
    else
        xc, vc = RotationCurve
        spl = Spline1D(ustrip.(uLen, xc), ustrip.(uVel, vc); k, bc)
        v = spl(ustrip.(uLen, R)) * uVel
        vel = rotational_velocity.(pos.x, pos.y, v, rotational_ratio)

        # Cancel out shifting
        v0 = mean(vel)
        vel = vel .- v0
        uVel = getuVel(units)
        vel = uconvert.(uVel, vel)
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