# using AstroPlot
using Random: seed!

@testset "SphericalSystem" begin
    r_s = 2.0u"kpc"
    MaxRadius = 3*r_s
    NumSamples = 100

    r = collect(LinRange(0.001u"kpc", MaxRadius, 100))
    dr = r[2] - r[1]

    shell_mass_func = x->1.0e8u"Msun/kpc^3"/(x/r_s)/(1+x/r_s)^2 * 4π * x^2

    config = SphericalSystem(STAR, NumSamples, shell_mass_func)
    particles = generate(config; MaxRadius)
    @test length(particles) == NumSamples
    # plot_makie(particles)


    shell_mass = shell_mass_func.(r)

    config = SphericalSystem(STAR, NumSamples, Dict("r"=>r, "m"=>shell_mass))
    particles = generate(config; MaxRadius)
    @test length(particles) == NumSamples
    # plot_makie(particles)
end

@testset "SphericalSystem: RotationCurve" begin
    r_s = 2.0u"kpc"
    shell_mass_func = x -> 1.0e8u"Msun/kpc^3" / (x/r_s) / (1 + x/r_s)^2 * 4π * x^2
    config = SphericalSystem(STAR, 30, shell_mass_func)

    rc_r = [0.5, 1.0, 2.0, 4.0, 6.0] .* u"kpc"
    rc_v = [50.0, 100.0, 150.0, 200.0, 220.0] .* u"km/s"
    particles = generate(config; MaxRadius = 6.0u"kpc", RotationCurve = (rc_r, rc_v))

    @test length(particles) == 30
    @test all(p -> !iszero(sqrt(p.Vel.x^2 + p.Vel.y^2)), particles)
end

@testset "SphericalSystem: VelocityDispersion" begin
    seed!(9)
    config = SphericalSystem(STAR, 200, x -> 1.0u"Msun/kpc")
    particles = generate(config; MaxRadius = 5.0u"kpc",
                         VelocityDispersion = 50.0u"km/s")

    vx = [ustrip(u"km/s", p.Vel.x) for p in particles]
    @test any(!iszero, vx)
    @test std(vx) ≈ 50.0 rtol = 0.25
    @test ustrip(u"km/s", norm(mean(particles.Vel))) < 10.0

    # Radial profile form
    rd = r -> 50.0u"km/s" * (r / 1.0u"kpc")
    p2 = generate(config; MaxRadius = 5.0u"kpc", VelocityDispersion = rd)
    @test all(!iszero(norm(p.Vel)) for p in p2)
end

@testset "SphericalSystem: collection propagation" begin
    config = SphericalSystem(GAS, 10, x -> 1.0u"Msun/kpc")
    particles = generate(config; MaxRadius = 5.0u"kpc")
    @test all(==(GAS), particles.Collection)
end

@testset "SphericalSystem: Base.show" begin
    config = SphericalSystem(STAR, 100, Dict("r"=>[1.0u"kpc"], "m"=>[1.0e8u"Msun/kpc"]))
    io = IOBuffer()
    show(io, config)
    @test !isempty(String(take!(io)))
end

@testset "SphericalSystem: error paths" begin
    r_s = 2.0u"kpc"
    MaxRadius = 3 * r_s
    NumSamples = 50
    shell_mass_func = x -> 1.0e8u"Msun/kpc^3" / (x/r_s) / (1 + x/r_s)^2 * 4π * x^2

    # MaxRadius = 0 → error
    @test_throws ErrorException generate(SphericalSystem(STAR, NumSamples, shell_mass_func); MaxRadius = 0.0u"kpc")

    # Dict with missing keys → error
    bad_dict = Dict("a" => 1, "b" => 2)
    @test_throws ErrorException generate(SphericalSystem(STAR, NumSamples, bad_dict); MaxRadius = MaxRadius)

    # DataFrame with missing properties → error
    bad_df = DataFrame(x = [1, 2, 3])
    @test_throws ErrorException generate(SphericalSystem(STAR, NumSamples, bad_df); MaxRadius = MaxRadius)

    # Unsupported type → error
    @test_throws ErrorException generate(SphericalSystem(STAR, NumSamples, "not a valid type"); MaxRadius = MaxRadius)
end

@testset "SphericalSystem: discrete mass_shell forms" begin
    r_s = 2.0u"kpc"
    MaxRadius = 3 * r_s
    NumSamples = 50
    r = collect(LinRange(0.001u"kpc", MaxRadius, 50))
    shell_mass = 1.0e8u"Msun/kpc^3" ./ (r ./ r_s) ./ (1 .+ r ./ r_s).^2 .* 4π .* r.^2

    # DataFrame with :r / :m columns.
    df = DataFrame(r = r, m = shell_mass)
    p_df = generate(SphericalSystem(STAR, NumSamples, df); MaxRadius)
    @test length(p_df) == NumSamples

    # Dict with either symbol or string keys.
    p_sym = generate(
        SphericalSystem(STAR, NumSamples, Dict(:r => r, :m => shell_mass));
        MaxRadius,
    )
    p_str = generate(
        SphericalSystem(STAR, NumSamples, Dict("r" => r, "m" => shell_mass));
        MaxRadius,
    )
    @test length(p_sym) == NumSamples
    @test length(p_str) == NumSamples

    # The table is interpreted in its own units, independent of the output unit.
    seed!(123)
    p_astro = generate(SphericalSystem(STAR, NumSamples, df), uAstro; MaxRadius)
    seed!(123)
    p_si = generate(SphericalSystem(STAR, NumSamples, df), uSI; MaxRadius)
    for (a, b) in zip(p_astro, p_si)
        @test ustrip(u"kpc", norm(a.Pos)) ≈ ustrip(u"kpc", norm(b.Pos))
    end
end
