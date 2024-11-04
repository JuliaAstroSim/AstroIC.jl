"""
$(TYPEDSIGNATURES)

Use `Optim.jl` to find the minimum value of `target` function.
"""
function minimum_func(target::Function, ic)
    result = Optim.optimize(target, ic)
    maximum_value = Optim.minimum(result)
    optimal_xy = Optim.minimizer(result)
    return maximum_value, optimal_xy
end

"""
$(TYPEDSIGNATURES)

Ref: https://en.wikipedia.org/wiki/Exponential_distribution
"""
function exponential_pdf(beta, x)
    if x > zero(x)
        y = exp(-x/beta)/beta
    else
        y = zero(x)
    end
    return y
end

"""
$(TYPEDSIGNATURES)

Ref: https://en.wikipedia.org/wiki/Exponential_distribution
"""
function exponential_cdf(beta, x)
    if x > zero(x)
        y = 1-exp(-x/beta)
    else
        y = zero(x)
    end
    return y
end

"""
$(TYPEDSIGNATURES)

Ref: https://en.wikipedia.org/wiki/Exponential_distribution
"""
function exponential_cdf_inv(beta, y)
    x = -beta*log(1-y)
end

"""
$(TYPEDSIGNATURES)

"""
function exponential_pdf_hole(beta, Rm, x)
    if x > zero(x)
        y = exp(-x/beta-Rm/x)/beta
    else
        y = zero(x)
    end
    return y
end

"""
$(TYPEDSIGNATURES)

"""
function sech2_pdf(x)
    return sech(x)^2
end

"""
$(TYPEDSIGNATURES)

"""
function sech2_cdf(x)
    return 0.5*tanh(x)+0.5
end

"""
$(TYPEDSIGNATURES)

"""
function sech2_cdf_inv(u)
    return atanh(2*u-1)
end