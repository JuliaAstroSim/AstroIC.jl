# using AstroPlot

@testset "SphericalSystem" begin
    r_s = 2.0u"kpc"
    MaxRadius = 10*r_s
    NumSamples = 1000

    r = collect(LinRange(0.001u"kpc", MaxRadius, 1000))
    dr = r[2] - r[1]

    shell_mass_func = x->1.0e8u"Msun/kpc^3"/(x/r_s)/(1+x/r_s)^2 * 4π * x^2
    
    config = SphericalSystem(STAR, NumSamples, shell_mass_func)
    particles = generate(config; MaxRadius)
    # plot_makie(particles)


    shell_mass = shell_mass_func.(r)

    config = SphericalSystem(STAR, NumSamples, Dict("r"=>r, "m"=>shell_mass))
    particles = generate(config; MaxRadius)
    # plot_makie(particles)
end