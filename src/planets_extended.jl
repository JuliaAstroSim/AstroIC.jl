# Extended Solar System generators
#
# Adds three functions on top of the existing `solarsystem(date)`:
#
#   1. `moon(date; body="Moon", parent="Earth")` — single body via JPL Horizons REST
#   2. `solar_system_full(date; include_moons=true)` — Sun + 8 planets + Pluto
#      + Moon + Galilean moons + Titan, all in one StructArray
#   3. `horizons_state_vector(body, date; center="500@10", step="1 d")` — raw
#      call to JPL Horizons REST API, returns (x, y, z, vx, vy, vz) in km & km/s
#
# Why call Horizons instead of using AstroLib?
#   - AstroLib only supports the 8 planets + Sun (no Moon position, no Galilean moons,
#     no Pluto, no spacecraft).
#   - For the Outreach blog series we want to show the Moon + Galilean moons + Titan
#     visually alongside the planets; the only way is to query JPL directly.
#
# Caching: `solar_system_full` stores results in a small in-memory cache keyed by
# the body's horizons_id + the date rounded to the nearest hour. Repeated calls
# during one session are free.

using Dates

# ---------------------------------------------------------------------------
# In-memory cache (per-process). Keyed by (body, center, step, exact epoch).
# Including center/step is essential: the same body has different state
# vectors in different reference frames.
# ---------------------------------------------------------------------------

const _HORIZONS_CACHE = Dict{Tuple{String,String,String,DateTime}, NTuple{6,Float64}}()

# ---------------------------------------------------------------------------
# Low-level: query JPL Horizons REST API
# ---------------------------------------------------------------------------

"""
    horizons_state_vector(body::AbstractString, date::DateTime;
                          center::AbstractString = "500@10",
                          step::AbstractString = "1 d") -> NTuple{6,Float64}

Query NASA JPL Horizons for the state vector of `body` at `date`, returning
`(x, y, z, vx, vy, vz)` in **km** and **km/s** in the **J2000 ecliptic frame**.

`body` can be any Horizons COMMAND id recognized by the system — for example
`"10"` (Sun), `"199"` (Mercury), `"301"` (Moon), `"501"` (Io), `"999"` (Pluto).

`center` is the coordinate origin:
- `"500@10"`  — Sun body center (heliocentric)
- `"500@399"` — Earth body center
- `"500@599"` — Jupiter body center
- `"500@699"` — Saturn body center
- `"500@0"`   — solar-system barycenter

`step` is the ephemeris step (defaults to `1 d`, so the call is cheap).

Results are cached per exact epoch, center and step (`_HORIZONS_CACHE`), so
repeated identical calls during a session are free.

## Example

```julia
julia> using AstroIC, Dates
julia> sv = AstroIC.horizons_state_vector("301", DateTime(2026, 7, 25, 0, 0); center="500@399")
# (geocentric position + velocity of the Moon in km, km/s)
```
"""
function horizons_state_vector(body::AbstractString, date::DateTime;
                               center::AbstractString = "500@10",
                               step::AbstractString = "1 d")
    cache_key = (String(body), String(center), String(step), date)
    if haskey(_HORIZONS_CACHE, cache_key)
        return _HORIZONS_CACHE[cache_key]
    end

    date_str = replace(Dates.format(date, "yyyy-mm-dd HH:MM"), " " => "%20")
    step_str = replace(step, " " => "%20")

    url = "https://ssd.jpl.nasa.gov/api/horizons.api?" *
          "format=text&" *
          "COMMAND='" * body * "'&" *
          "OBJ_DATA=NO&MAKE_EPHEM=YES&" *
          "EPHEM_TYPE=VECTORS&" *
          "CENTER='" * center * "'&" *
          "START_TIME='" * date_str * "'&" *
          "STOP_TIME='" * date_str * "'&" *
          "STEP_SIZE='" * step_str * "'&" *
          "VEC_TABLE=2&" *
          "REF_PLANE='ECLIPTIC'&" *
          "OUT_UNITS='KM-S'"

    txt = _http_get(url)

    sv = _parse_horizons_state_vector(txt, body, date)
    _HORIZONS_CACHE[cache_key] = sv
    return sv
end

# Simple HTTP GET helper using the Downloads stdlib. The response is written
# into an in-memory buffer; note that `Downloads.download(url)` alone returns a
# temporary file path, not the body text.
function _http_get(url::AbstractString)
    io = IOBuffer()
    try
        Downloads.download(url, io; timeout = 30)
        return String(take!(io))
    catch e
        error("""
        Could not reach JPL Horizons at $url.
        Reason: $(sprint(showerror, e))

        If you are offline, this is expected. Re-run when you have network access,
        or pre-populate `_HORIZONS_CACHE` with cached values from a previous run.
        """)
    end
end

# Parse the plain-text Horizons response. The block we care about looks like:
#
#  $$SOE
#   2460580.500000000 = A.D. 2026-Jul-25 00:00:00.0000  TDB
#        X = 1.234567890123456E+05 Y = ... Z = ...
#        VX= ... VY= ... VZ= ...
#  $$EOE
#
# We extract the first X/Y/Z and VX/VY/VZ pair after `$$SOE`.

function _parse_horizons_state_vector(txt::AbstractString, body::AbstractString, date::DateTime)
    # Find the SOE/EOE markers
    soe_idx = findfirst("SOE", txt)
    eoe_idx = findfirst("EOE", txt)
    if soe_idx === nothing || eoe_idx === nothing
        error("Could not parse Horizons response for body=$(body), date=$(date). " *
              "Raw response (first 500 chars):\n$(first(txt, 500))")
    end

    # Extract the block between markers
    block = txt[nextind(txt, last(soe_idx)):prevind(txt, first(eoe_idx))]

    # Extract X, Y, Z, VX, VY, VZ
    x = _extract_number(block, "X =")
    y = _extract_number(block, "Y =")
    z = _extract_number(block, "Z =")
    vx = _extract_number(block, "VX=")
    vy = _extract_number(block, "VY=")
    vz = _extract_number(block, "VZ=")

    return (x, y, z, vx, vy, vz)
end

function _extract_number(block::AbstractString, key::AbstractString)
    idx = findfirst(key, block)
    if idx === nothing
        error("Could not find key '$(key)' in Horizons response block:\n$(first(block, 200))")
    end
    # Read the first token after the key (skipping whitespace)
    rest = block[nextind(block, last(idx)):end]
    m = match(r"^\s*([-+]?\d+\.?\d*(?:[Ee][-+]?\d+)?)", rest)
    if m === nothing
        error("Could not parse number for key '$(key)' in: $(first(rest, 60))")
    end
    return parse(Float64, m.captures[1])
end

# ---------------------------------------------------------------------------
# Single-body wrappers
# ---------------------------------------------------------------------------

"""
    moon(date::DateTime; body::AbstractString = "Moon", parent::AbstractString = "Earth")

Return a `Star` (PhysicalParticle) for `body` (default `"Moon"`) at `date`,
positioned relative to its `parent` body.

Currently supported `body` values (case-sensitive, see [`PLANET_TABLE`](@ref)):
`"Moon"`, `"Io"`, `"Europa"`, `"Ganymede"`, `"Callisto"`, `"Titan"`.

By default the parent is taken from [`PLANET_TABLE`](@ref) (`Earth` for the
Moon, `Jupiter` for the Galilean moons, `Saturn` for Titan), but it can be
overridden with the `parent` keyword. The returned `Star` uses SI units.
"""
function moon(date::DateTime; body::AbstractString = "Moon",
              parent::Union{Nothing,AbstractString} = nothing)
    info = planet_info(body)
    parent_name = isnothing(parent) ? info.parent : parent
    isempty(parent_name) && error("Body `$(body)` has no parent; specify `parent`.")
    parent_center = "500@" * planet_horizons_id(parent_name)

    # Step of 1 hour is fine — the moon moves only ~0.04° in an hour
    x, y, z, vx, vy, vz = horizons_state_vector(info.horizons_id, date;
                                                 center = parent_center,
                                                 step   = "1 h")

    p = Star(uSI, id = -1)  # id = -1 means "moon" (planet ids are 1..8, sun is 0)
    p = setproperty!!(p, :Pos, PVector(
        uconvert(u"m", u"km" * x),
        uconvert(u"m", u"km" * y),
        uconvert(u"m", u"km" * z),
    ))
    p = setproperty!!(p, :Vel, PVector(
        uconvert(u"m/s", u"km/s" * vx),
        uconvert(u"m/s", u"km/s" * vy),
        uconvert(u"m/s", u"km/s" * vz),
    ))
    p = setproperty!!(p, :Mass, info.M_kg * u"kg")
    return p
end

# ---------------------------------------------------------------------------
# Full Solar System (sun + 8 planets + Pluto + selected moons)
# ---------------------------------------------------------------------------

"""
    solar_system_full(date::DateTime;
                       include_moons::Bool = true,
                       include_pluto::Bool = true) -> StructArray

Like [`solarsystem`](@ref), but **complete**: sun + 8 planets + Pluto + Moon +
Galilean moons (Io, Europa, Ganymede, Callisto) + Titan — 16 bodies in total.

Each body is a `Star` (PhysicalParticle) with fields:
- `Pos`     : heliocentric Cartesian position (SI units, m)
- `Vel`     : heliocentric velocity (SI units, m/s)
- `Mass`    : SI kg

Bodies are in **heliocentric ecliptic J2000 Cartesian** coordinates, even the
moons. Moon position is computed as `parent_position + parent-centered offset`
so all bodies share the same frame.

> **Performance**: the first call performs several sequential Horizons
> requests. Results are cached per-process (keyed by body, center, step and
> exact epoch), so subsequent identical calls are free.

```julia
julia> using AstroIC, Dates
julia> bodies = solar_system_full(DateTime(2026, 7, 25))
```
"""
function solar_system_full(date::DateTime;
                           include_moons::Bool = true,
                           include_pluto::Bool = true)
    # 1) Start with the existing 8 planets + Sun
    bodies = solarsystem(date)

    if include_pluto
        push!(bodies, _planet_particle(date, "Pluto"))
    end

    if include_moons
        # Every body with a non-Sun parent in PLANET_TABLE is a moon. Its
        # state vector is returned RELATIVE to that parent, so we add the
        # parent's heliocentric position to convert to heliocentric.
        for info in PLANET_TABLE
            (isempty(info.parent) || info.parent == "Sun") && continue
            try
                push!(bodies, _moon_particle_heliocentric(date, info.name, info.parent))
            catch e
                @warn "Failed to fetch $(info.name): $(e)"
            end
        end
    end

    return bodies
end

# Helper: build a heliocentric PhysicalParticle for a body whose state is
# returned relative to a parent.
function _moon_particle_heliocentric(date::DateTime, body::AbstractString, parent::AbstractString)
    # Heliocentric position = parent_position (heliocentric) + offset (parent-centered)
    # Heliocentric velocity  = parent_velocity  + offset_velocity

    # 1) Parent state in heliocentric coordinates.
    parent_id = planet_horizons_id(parent)
    parent_x, parent_y, parent_z, parent_vx, parent_vy, parent_vz =
        horizons_state_vector(parent_id, date; center = "500@10", step = "1 h")

    # 2) Offset (relative to parent), in km & km/s
    rel_x, rel_y, rel_z, rel_vx, rel_vy, rel_vz =
        horizons_state_vector(planet_horizons_id(body), date;
                              center = "500@" * parent_id,
                              step   = "1 h")

    x_km  = parent_x + rel_x
    y_km  = parent_y + rel_y
    z_km  = parent_z + rel_z
    vx_kms = parent_vx + rel_vx
    vy_kms = parent_vy + rel_vy
    vz_kms = parent_vz + rel_vz

    p = Star(uSI, id = -1)
    p = setproperty!!(p, :Pos, PVector(u"m" * (x_km  * 1000),
                                       u"m" * (y_km  * 1000),
                                       u"m" * (z_km  * 1000)))
    p = setproperty!!(p, :Vel, PVector(u"m/s" * (vx_kms * 1000),
                                       u"m/s" * (vy_kms * 1000),
                                       u"m/s" * (vz_kms * 1000)))
    p = setproperty!!(p, :Mass, planet_info(body).M_kg * u"kg")
    return p
end

# Helper: build a heliocentric particle for Pluto (no parent — direct query)
function _planet_particle(date::DateTime, name::AbstractString)
    info = planet_info(name)
    x, y, z, vx, vy, vz = horizons_state_vector(info.horizons_id, date;
                                                  center = "500@10",
                                                  step   = "1 h")
    p = Star(uSI, id = -1)
    p = setproperty!!(p, :Pos, PVector(u"m" * (x  * 1000),
                                       u"m" * (y  * 1000),
                                       u"m" * (z  * 1000)))
    p = setproperty!!(p, :Vel, PVector(u"m/s" * (vx * 1000),
                                       u"m/s" * (vy * 1000),
                                       u"m/s" * (vz * 1000)))
    p = setproperty!!(p, :Mass, info.M_kg * u"kg")
    return p
end