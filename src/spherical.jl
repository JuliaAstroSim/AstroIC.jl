"""
$(TYPEDEF)
$(TYPEDFIELDS)

`mass_shell` support two types:
- `Dict` or `DataFrames`: discrete radii `:r` and shell masses `:m`
- `function`: a radial function of shell mass
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

"""
function generate(config::SphericalSystem, units = uAstro;
    MaxRadius = 0.0u"kpc",
)
    uLen = getuLength(units)
    uMass = getuMass(units)
    uVel = getuVel(units)
    NumSamples = config.NumSamples

    if iszero(MaxRadius)
        error("Keyword `MaxRadius` should not be zero")
    end

    if config.mass_shell isa Function
        R = rejection_sampling(config.mass_shell, MaxRadius, NumSamples)
    elseif config.mass_shell isa Dict
        if haskey(config.mass_shell, "r") && haskey(config.mass_shell, "m")
            R = rejection_sampling(config.mass_shell["r"], config.mass_shell["m"], MaxRadius, NumSamples)
        else
            error("`mass_shell` of `SphericalSystem` must be `Dict` or `DataFrame` (with fields `:r`, `:m`) or `Function`")
        end
    elseif config.mass_shell isa DataFrame
        if hasproperty(config.mass_shell, :r) && hasproperty(config.mass_shell, :m)
            R = rejection_sampling(config.mass_shell.r, config.mass_shell.m, MaxRadius, NumSamples)
        else
            error("`mass_shell` of `SphericalSystem` must be `Dict` or `DataFrame` (with fields `:r`, `:m`) or `Function`")
        end
    else
        error("`mass_shell` of `SphericalSystem` must be `Dict` or `DataFrame` (with fields `:r`, `:m`) or `Function`")
    end

    # Compute total mass
    x = collect(LinRange(MaxRadius/10000, MaxRadius, 10000))
    dx = x[2] - x[1]
    if config.mass_shell isa Function
        shell_masses = config.mass_shell.(x)
        TotalMass = PhysicalParticles.NumericalIntegration.integrate(x, shell_masses)
    elseif config.mass_shell isa Dict || config.mass_shell isa DataFrame
        spl = Spline1D(ustrip.(uLen, config.mass_shell["r"]), ustrip.(uMass/uLen, config.mass_shell["m"]))
        TotalMass = PhysicalParticles.NumericalIntegration.integrate(x, spl(ustrip.(uLen, x))*uMass/uLen)
    end

    pos = rand_pos_3d.(R)

    #TODO: vel
    vel = [PVector(uVel) for i in 1:NumSamples]

    # Packing
    particles = StructArray(Star(units, id = i) for i in 1:NumSamples)
    assign_particles(particles, :Pos, pos)
    assign_particles(particles, :Vel, vel)

    Mmean = TotalMass / NumSamples
    assign_particles(particles, :Mass, Mmean)

    return particles
end