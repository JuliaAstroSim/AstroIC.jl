getuVel(::Nothing) = nothing
getuVel(units::Array) = units[1] / units[2]

struct PlummerStarCluster{T<:AbstractParticleCollection} <: InitialConditionConfig{T}
    collection::T
    NumSamples::Integer

    VirialRadius::Number

    TotalMass::Number
    G::Number

    model::GravModel
end

function PlummerStarCluster(;
        collection = STAR(),
        NumSamples = 1000,
        VirialRadius = 0.010u"kpc",
        TotalMass = 1.0e5u"Msun",
        G = 4.498502151469553e-6u"kpc^3 / Msun / Gyr^2",

        model = Newton(),
    )
    
    return PlummerStarCluster(
        collection,
        NumSamples,
        VirialRadius,
        TotalMass,
        G,
        model
    )
end

function Base.show(io::IO, config::PlummerStarCluster)
    print(
        io,
        "Config of Plummer Star Cluster Initial Conditions:",
        "\n          Gravity Model: ", typeof(config.model),
        "\n    Particle collection: ", config.collection,
        "\n      Number of Samples: ", config.NumSamples,
        "\n          Virial radius: ", config.VirialRadius,
        "\n             Total Mass: ", config.TotalMass,
        "\n       Gravity constant: ", config.G,
    )
end

@inline plummer_pdfr_inv(VirialRadius::Number, rnd::AbstractFloat) = VirialRadius/(rnd^(-2.0/3.0)-1.0)^0.5

function rand_plummerpos(r::Number, VirialRadius::Number)
    eta, gamma = rand(2)

    x = r * sin(acos(2.0eta-1.0)) * cos(2.0pi*gamma)
    y = r * sin(acos(2.0eta-1.0)) * sin(2.0pi*gamma)
    z = r * (2.0eta-1.0)
    return PVector(x, y, z)
end

function rand_plummerpos(r::Array{T,N}, VirialRadius::Number) where T<:Number where N
    return [rand_plummerpos(i, VirialRadius) for i in r]
end

function plummer_vel_sigma2(r::Number, VirialRadius::Number, Mass::Number, G::Number, ::Newton)
    return G * Mass / sqrt(r^2 + VirialRadius^2) / 6.0
end

function plummer_vel_sigma2(r::Number, VirialRadius::Number, Mass::Number, G::Number, ::MOND)
    return G * Mass / sqrt(r^2 + VirialRadius^2) / 2.0
end

function rand_plummervel(r::Number, VirialRadius::Number, Mass::Number, G::Number, model::GravModel)
    v = sqrt(plummer_vel_sigma2(r, VirialRadius, Mass, G, model))
    return PVector(randn() * v, randn() * v, randn() * v)
end

function rand_plummervel(r::Array{T,N}, VirialRadius::Number, Mass::Number, G::Number, model::GravModel) where T<:Number where N
    return [rand_plummervel(i, VirialRadius, Mass, G, model) for i in r]
end

function generate(config::PlummerStarCluster, units = uAstro;
                  MaxRadius = 5 * config.VirialRadius,
                  )
    println(config)

    uLength = getuLength(units)

    NumSamples = config.NumSamples
    VirialRadius = uconvert(uLength, config.VirialRadius)

    # generate radii
    r = plummer_pdfr_inv.(VirialRadius, rand(NumSamples))
    
    # MaxRadius
    if MaxRadius > VirialRadius
        for i in eachindex(r)
            while r[i] > MaxRadius
                r[i] = plummer_pdfr_inv(VirialRadius, rand())
            end
        end
    end

    # Sampling
    pos = rand_plummerpos(r, VirialRadius)
    vel = rand_plummervel(r, VirialRadius, config.TotalMass, config.G, config.model)

    # Cancel out shifting
    v0 = mean(vel)
    vel = vel .- v0
    uVel = getuVel(units)
    vel = uconvert.(uVel, vel)
    
    # Packing
    particles = [Star(uAstro, id = i) for i in 1:NumSamples]
    assign_particles(particles, :Pos, pos)
    assign_particles(particles, :Vel, vel)
    
    Mmean = config.TotalMass / NumSamples
    for p in particles
        p.Mass = Mmean
    end

    return particles
end