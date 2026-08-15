# Planetary / dwarf-planet constants
#
# Reference values for the 8 planets, Pluto, the Moon, and the major moons that
# the Outreach blog series highlights (Io, Europa, Ganymede, Callisto, Titan).
#
# Sources:
#   - IAU 2015 Resolution B3 (nominal values)
#   - NASA Planetary Fact Sheet (https://nssdc.gsfc.nasa.gov/planetary/factsheet/)
#   - NASA JPL Horizons (https://ssd.jpl.nasa.gov/horizons/) — Horizons COMMAND id
#
# Conventions:
#   - a_AU         : semi-major axis in AU
#   - e            : orbital eccentricity (dimensionless)
#   - i_deg        : orbital inclination to ecliptic (degrees)
#   - P_year       : orbital period (years)
#   - r_eq_km      : mean equatorial radius (km)
#   - M_kg         : mass (kg)
#   - albedo       : geometric albedo (dimensionless)
#   - color        : (R, G, B) ∈ [0,1] from NASA style guide (color-blind friendly)
#   - horizons_id  : NASA JPL Horizons COMMAND id (string) for `horizons_state(...)`
#   - parent       : name of the body this object orbits (empty string for heliocentric)

const PLANET_TABLE = [
    (name = "Sun",       horizons_id = "10",  parent = "",        a_AU = 0.0,      e = 0.0,    i_deg = 0.0,    P_year = 0.0,       r_eq_km = 695700.0, M_kg = 1.989e30,  albedo = 1.00, color = (1.000, 0.965, 0.840)),
    (name = "Mercury",   horizons_id = "199", parent = "Sun",     a_AU = 0.38710,  e = 0.2056, i_deg = 7.005,  P_year = 0.24085,   r_eq_km = 2439.7,   M_kg = 3.3011e23, albedo = 0.14, color = (0.660, 0.660, 0.660)),
    (name = "Venus",     horizons_id = "299", parent = "Sun",     a_AU = 0.72333,  e = 0.0068, i_deg = 3.395,  P_year = 0.61520,   r_eq_km = 6051.8,   M_kg = 4.8675e24, albedo = 0.72, color = (0.972, 0.870, 0.732)),
    (name = "Earth",     horizons_id = "399", parent = "Sun",     a_AU = 1.00000,  e = 0.0167, i_deg = 0.000,  P_year = 1.00000,   r_eq_km = 6378.1,   M_kg = 5.972e24,  albedo = 0.30, color = (0.255, 0.412, 0.918)),
    (name = "Moon",      horizons_id = "301", parent = "Earth",   a_AU = 0.00257,  e = 0.0549, i_deg = 5.145,  P_year = 0.07480,   r_eq_km = 1737.4,   M_kg = 7.342e22,  albedo = 0.12, color = (0.831, 0.812, 0.769)),
    (name = "Mars",      horizons_id = "499", parent = "Sun",     a_AU = 1.52366,  e = 0.0934, i_deg = 1.851,  P_year = 1.88090,   r_eq_km = 3389.5,   M_kg = 6.417e23,  albedo = 0.25, color = (0.756, 0.323, 0.219)),
    (name = "Jupiter",   horizons_id = "599", parent = "Sun",     a_AU = 5.20336,  e = 0.0484, i_deg = 1.305,  P_year = 11.86200,  r_eq_km = 69911.0,  M_kg = 1.898e27,  albedo = 0.52, color = (0.847, 0.706, 0.533)),
    (name = "Io",        horizons_id = "501", parent = "Jupiter", a_AU = 0.00282,  e = 0.0041, i_deg = 0.036,  P_year = 0.00485,   r_eq_km = 1821.6,   M_kg = 8.932e22,  albedo = 0.63, color = (0.972, 0.870, 0.624)),
    (name = "Europa",    horizons_id = "502", parent = "Jupiter", a_AU = 0.00448,  e = 0.0094, i_deg = 0.469,  P_year = 0.00972,   r_eq_km = 1560.8,   M_kg = 4.799e22,  albedo = 0.67, color = (0.870, 0.847, 0.780)),
    (name = "Ganymede",  horizons_id = "503", parent = "Jupiter", a_AU = 0.00716,  e = 0.0011, i_deg = 0.177,  P_year = 0.01959,   r_eq_km = 2634.1,   M_kg = 1.482e23,  albedo = 0.43, color = (0.660, 0.624, 0.580)),
    (name = "Callisto",  horizons_id = "504", parent = "Jupiter", a_AU = 0.01259,  e = 0.0074, i_deg = 0.192,  P_year = 0.04569,   r_eq_km = 2410.3,   M_kg = 1.076e23,  albedo = 0.22, color = (0.435, 0.412, 0.385)),
    (name = "Saturn",    horizons_id = "699", parent = "Sun",     a_AU = 9.53707,  e = 0.0542, i_deg = 2.484,  P_year = 29.45710,  r_eq_km = 58232.0,  M_kg = 5.683e26,  albedo = 0.47, color = (0.918, 0.831, 0.658)),
    (name = "Titan",     horizons_id = "606", parent = "Saturn",  a_AU = 0.00817,  e = 0.0288, i_deg = 0.349,  P_year = 0.04359,   r_eq_km = 2574.7,   M_kg = 1.345e23,  albedo = 0.22, color = (0.870, 0.624, 0.219)),
    (name = "Uranus",    horizons_id = "799", parent = "Sun",     a_AU = 19.19126, e = 0.0472, i_deg = 0.770,  P_year = 84.01100,  r_eq_km = 25362.0,  M_kg = 8.681e25,  albedo = 0.51, color = (0.706, 0.870, 0.918)),
    (name = "Neptune",   horizons_id = "899", parent = "Sun",     a_AU = 30.06896, e = 0.0086, i_deg = 1.770,  P_year = 164.79000, r_eq_km = 24622.0,  M_kg = 1.024e26,  albedo = 0.41, color = (0.412, 0.580, 0.918)),
    (name = "Pluto",     horizons_id = "999", parent = "Sun",     a_AU = 39.48169, e = 0.2488, i_deg = 17.160, P_year = 247.94000, r_eq_km = 1188.3,   M_kg = 1.303e22,  albedo = 0.30, color = (0.847, 0.706, 0.580)),
]

"""
    planet_info(name::AbstractString)

Return the entry of [`PLANET_TABLE`](@ref) whose `name` matches `name` (case-sensitive).
Throws `ArgumentError` if no match.

```julia
julia> planet_info("Jupiter")
(name = "Jupiter", horizons_id = "599", parent = "Sun", a_AU = 5.20336, ...)
```
"""
function planet_info(name::AbstractString)
    for p in PLANET_TABLE
        if p.name == name
            return p
        end
    end
    throw(ArgumentError("Unknown planet: $name. Available: $([p.name for p in PLANET_TABLE])"))
end

"""
    planet_names() -> Vector{String}

Return the names of every body in [`PLANET_TABLE`](@ref).
"""
planet_names() = [p.name for p in PLANET_TABLE]

"""
    planet_color(name::AbstractString) -> Tuple{Float64,Float64,Float64}

Return the (R, G, B) color of the body `name`. Used for plotting.
"""
function planet_color(name::AbstractString)
    return planet_info(name).color
end

"""
    planet_horizons_id(name::AbstractString) -> String

Return the NASA JPL Horizons COMMAND id (string) of the body `name`.
For example: `planet_horizons_id("Moon") == "301"`, `planet_horizons_id("Jupiter") == "599"`.
"""
planet_horizons_id(name::AbstractString) = planet_info(name).horizons_id