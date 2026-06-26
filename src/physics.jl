"""
    vmean(T::Number, mass::Number, units = uAstro)

Thermal (Maxwell–Boltzmann) mean speed ``\\sqrt{8 k_B T / (\\pi m)}`` of a particle
of mass `mass` at temperature `T`, converted to the velocity unit of `units`.

# Arguments
- `T`: temperature (with units, e.g. `300u"K"`)
- `mass`: particle mass (with units, e.g. `Constant().m_p` or `1.0u"Msun"`)
- `units`: target unit system (default `uAstro`); must support `getuVel`

# Returns
A velocity in the velocity unit of `units` (e.g. `kpc/Gyr` for `uAstro`,
`m/s` for `uSI`).

# Example
```julia
v = vmean(300u"K", Constant().m_p, uSI)   # ≈ 1.15 km/s for a proton
```
"""
function vmean(T::Number, mass::Number, units = uAstro)
    return uconvert(getuVel(units), sqrt(8 * Constant(units).k_B * T / pi / mass))
end

"""
    vmean2(T::Number, mass::Number, units = uAstro)

Squared thermal mean speed ``8 k_B T / (\\pi m)`` of a particle of mass `mass`
at temperature `T`, converted to the squared-velocity unit of `units`.

See also [`vmean`](@ref).
"""
function vmean2(T::Number, mass::Number, units = uAstro)
    return uconvert(getuVel(units), 8 * Constant(units).k_B * T / pi / mass)
end