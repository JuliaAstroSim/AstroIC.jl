function vmean(T::Number, mass::Number, units = uAstro)
    return uconvert(getuVel(units), sqrt(8 * Constant(units).k_B * T / pi / mass))
end

function vmean2(T::Number, mass::Number, units = uAstro)
    return uconvert(getuVel(units), 8 * Constant(units).k_B * T / pi / mass)
end



"""
function softlength

    return recommended softening length
"""
function softlength(V::Number, N::Int64)
    return cbrt(V / N)
end