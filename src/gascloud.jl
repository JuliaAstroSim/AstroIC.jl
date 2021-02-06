mutable struct GasCloud{I, Len, Den, Tmp, MASS} <: InitialConditionConfig
    collection::Collection

    Radius::Len
    rho0::Den
    T::Tmp
    ParticleMass::MASS

    Nx::I
    Ny::I
    Nz::I
end

"""
GasCloud

## fields

- `collection`    particle type
- `Radius`        Radius of cloud
- `rho0`          Central Density
- `T`             Gas Temperature
- `Nx, Ny, Nz`    Grid Resolution
"""
function GasCloud(;
        collection::Collection = GAS,
        Radius::Number = 20u"kpc",
        rho0::Number = 1250u"Msun/kpc^3",
        T::Number = 300u"K",
        ParticleMass::Number = Constant().m_p,
        Nx::Int = 11,
        Ny::Int = 11,
        Nz::Int = 11,
    )
    return GasCloud(
        collection,
        Radius,
        rho0,
        T,
        ParticleMass,
        Nx, Ny, Nz,
    )
end

function Base.show(io::IO, config::GasCloud)
    print(io,
        "Config of Gad Cloud Initial Conditions:",
        "\n   Particle Collection: ", config.collection,
        "\n                Radius: ", config.Radius,
        "\n       Central Density: ", config.rho0,
        "\n       Gas Temperature: ", config.T,
        "\n       Grid Resolution: X = ", config.Nx, ", Y = ", config.Ny, ", Z = ", config.Nz,
    )
end

function gridpoints(R::Number, Nx::Int, Ny::Int, Nz::Int)
    x = zeros(Nx, Ny, Nz) * R
    y = zeros(Nx, Ny, Nz) * R
    z = zeros(Nx, Ny, Nz) * R

    Lx = 2.0 * R / (Nx - 1)
    Ly = 2.0 * R / (Ny - 1)
    Lz = 2.0 * R / (Nz - 1)

    rx = collect(-R:Lx:R)
    for j in 1:Nz
        for i in 1:Ny
            x[:,i,j] = rx
        end
    end

    ry = collect(-R:Ly:R)
    for j in 1:Nz
        for i in 1:Ny
            y[:,i,j] .= ry[i]
        end
    end

    rz = collect(-R:Lz:R)
    for j in 1:Nz
        z[:,:,j] .= rz[j]
    end

    return x, y, z
end

function generate(config::GasCloud, units = uAstro;
                  constants::Constant = Constant(units))
    println(config)

    R = config.Radius
    Nx = config.Nx
    Ny = config.Ny
    Nz = config.Nz
    Lx = 2.0 * R / (Nx - 1)
    Ly = 2.0 * R / (Ny - 1)
    Lz = 2.0 * R / (Nz - 1)

    x, y, z = gridpoints(config.Radius, config.Nx, config.Ny, config.Nz)
    pos = PVector.(x, y, z)

    v = vmean(config.T, config.ParticleMass, units)
    vrand = randn_pvector(Nx * Ny * Nz)
    vel = normalize.(vrand) * v
    
    data = Dict("gases" => empty([SPHGas(units)]))
    id = 1
    for i in 1:length(x)
        @inbounds r2 = pos[i] * pos[1]
        if iszero(r2)
            @inbounds r2 = pos[i+1] * pos[i+1]
        end

        if r2 <= R^2
            mass = config.rho0 * R^2 / r2 * Lx * Ly * Lz
            @inbounds push!(data["gases"], setproperties!!(SPHGas(units; id), Pos = pos[i], Vel = vel[i], Mass = mass))
            id += 1
        end
    end
    return data
end