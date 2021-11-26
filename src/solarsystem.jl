# Planets in Solar System
"""
    function helio2xyz(jd, num)

Convert heliocentric coordinates (in unit `AU`) to Cartesian coordinates in unit `m`

Returns a `PhysicalParticles::PVector`
"""
function helio2xyz(jd, num)
    r, theta, phi = helio(jd, num, true)
    return PVector(
        uconvert(u"m", r * u"AU" * sin(theta) * cos(phi)),
        uconvert(u"m", r * u"AU" * sin(theta) * sin(phi)),
        uconvert(u"m", r * u"AU" * cos(theta)),
    )
end

"""
    function solarsystem(date::DateTime)
    function solarsystem(date::Real)

Generate initial conditions of Solar System at desinated date
"""
function solarsystem(date::Real)
    # Planets
    #coords = [planet_coords(date, i) for i in 1:8]
    coords = [helio2xyz(date, i) for i in 1:8]
    
    names = [AstroLib.record[i] for i in 1:8]
    masses = [planets[names[i]].mass * u"kg" for i in 1:8]
    #axises = [planets[names[i]].axis * u"m" for i in 1:8]
    #periods = [planets[names[i]].period * u"s" for i in 1:8]
    #inclinations = [planets[names[i]].inc for i in 1:8]    # Degree
    #eccentricities = [planets[names[i]].ecc * u"s" for i in 1:8]

    coords1 = [helio2xyz(date + 1.0 / 86400.0, i) for i in 1:8] # Positions at next second
    vels = (coords1 - coords) / 1.0u"s"                         # Velocity is simply dx/dt

    particles = StructArray(Star(uSI, id = i) for i in 1:8)
    particles.Pos .= coords
    particles.Vel .= vels
    particles.Mass .= masses

    # Sun
    sun = Star(uSI)
    sun = setproperties!!(sun, Mass = uconvert(u"kg", 1.0u"Msun"))
    push!(particles, sun)
    #sunRA, sunDEC = sunpos(date)[1:2]
    #sunDistance = uconvert(u"m", 1.0u"AU")

    return particles
end

solarsystem(date::DateTime) = solarsystem(jdcnv(date))