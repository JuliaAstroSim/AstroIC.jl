"""
$(TYPEDSIGNATURES)

Use `Optim.jl` to find the minimum value of `target` function.
"""
function minimum_func(target::Function, ic)
    result = Optim.optimize(target, ic)
    # @show result
    minimum_value = Optim.minimum(result)
    optimal_xy = Optim.minimizer(result)
    return minimum_value, optimal_xy
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

"""
$(TYPEDSIGNATURES)

Sample from the discrete pdf using rejection method

!!! attentions
    - Data must be unitless, and have no NaN
    - `x` should start from 0
    - the maximum of `pdf` would be the width of sampling box
"""
function rejection_sampling(x, pdf, rMax, NumSamples)
    spl = Spline1D(x, pdf)
    pdf_maximum = maximum(pdf)

    R = eltype(pdf)[]
    sizehint!(R, NumSamples)

    while length(R) < NumSamples
        R_rand = rand() * rMax
        if rand() < spl(R_rand) / pdf_maximum
            push!(R, R_rand)
        end
    end
    
    return R
end