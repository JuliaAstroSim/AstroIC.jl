"""
$(TYPEDFIELDS)
"""
mutable struct PlummerStarCluster{I, Len, MASS, GM} <: InitialConditionConfig
    collection::Collection
    "amount of particles"
    NumSamples::I
    VirialRadius::Len
    "mass are seperated equally to all particles"
    TotalMass::MASS
    "gravity model [ Newton | MOND ]"
    model::GM
end

"""
$(TYPEDSIGNATURES)
"""
function PlummerStarCluster(;
        collection::Collection = STAR,
        NumSamples::Int64 = 1000,
        VirialRadius::Number = 0.010u"kpc",
        TotalMass::Number = 1.0e5u"Msun",

        model::GravityModel = Newton(),
    )
    
    return PlummerStarCluster(
        collection,
        NumSamples,
        VirialRadius,
        TotalMass,
        model
    )
end

function Base.show(io::IO, config::PlummerStarCluster)
    print(
        io,
        "Config of Plummer Star Cluster Initial Conditions:",
        "\n          Gravity Model: ", typeof(config.model),
        "\n    Particle Collection: ", config.collection,
        "\n      Number of Samples: ", config.NumSamples,
        "\n          Virial Radius: ", config.VirialRadius,
        "\n             Total Mass: ", config.TotalMass,
    )
end

@inline plummer_pdfr_inv(VirialRadius::Number, rnd::AbstractFloat) = VirialRadius/(rnd^(-2.0/3.0)-1.0)^0.5

function plummer_vel_sigma2(r::Number, VirialRadius::Number, Mass::Number, G::Number, ::Newton)
    return G * Mass / sqrt(r^2 + VirialRadius^2) / 6.0
end

function plummer_vel_sigma2(r::Number, VirialRadius::Number, Mass::Number, G::Number, ::MOND1983Milgrom)
    return G * Mass / sqrt(r^2 + VirialRadius^2) / 2.0
end

function rand_plummervel(r::Number, VirialRadius::Number, Mass::Number, G::Number, model::GravityModel)
    v = sqrt(plummer_vel_sigma2(r, VirialRadius, Mass, G, model))
    return PVector(randn() * v, randn() * v, randn() * v)
end

function rand_plummervel(r::Array{T,N}, VirialRadius::Number, Mass::Number, G::Number, model::GravityModel) where T<:Number where N
    return [rand_plummervel(i, VirialRadius, Mass, G, model) for i in r]
end

"""
$(TYPEDSIGNATURES)

Generate initial conditions of Plummer model

Automatically promote to `Measurement` if the parameters in config has `Measurement` type

# Keywords

- `MaxRadius`: resample particles outside the interested radius. Default is 5 * `VirialRadius`. Set to zero to avoid cutting off.

$_common_keywords
"""
function generate(config::PlummerStarCluster, units = uAstro;
                  constants::Constant = Constant(units),
                  MaxRadius = 5 * config.VirialRadius,
                  )
    uLength = getuLength(units)

    NumSamples = config.NumSamples
    VirialRadius = uconvert(uLength, config.VirialRadius)

    # generate radii
    if MaxRadius < VirialRadius
        @warn "`MaxRadius` is smaller than `VirialRadius`, this may cause unphyiscal errors!"
    end

    r = plummer_pdfr_inv.(VirialRadius, rand(NumSamples))
    for i in eachindex(r)
        while r[i] > MaxRadius
            r[i] = plummer_pdfr_inv(VirialRadius, rand())
        end
    end

    # Sampling
    pos = rand_pos_3d.(r)
    vel = rand_plummervel(r, VirialRadius, config.TotalMass, constants.G, config.model)

    # Cancel out shifting
    v0 = mean(vel)
    vel = vel .- v0
    uVel = getuVel(units)
    vel = uconvert.(uVel, vel)
    
    # Packing
    Promote = ismeasurement(config.VirialRadius) || ismeasurement(config.TotalMass)
    Promote && @info "Promoting to `Measurement`"
    particles = StructArray(Star(units, id = i, Measurement=Promote) for i in 1:NumSamples)
    assign_particles(particles, :Pos, pos)
    assign_particles(particles, :Vel, vel)
    
    Mmean = config.TotalMass / NumSamples
    assign_particles(particles, :Mass, Mmean)

    return particles
end