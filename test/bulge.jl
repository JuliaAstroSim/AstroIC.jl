@testset "Bulge" begin
    config = Bulge(;
        collection = STAR,
        NumSamples = 100,
        TotalMass = 8.57e9u"Msun",
        ScaleRadius = 0.075u"kpc",
        CutRadius = 2.1u"kpc",
        q = 0.5,
        α = 1.8,
    )

    data = generate(
        config,
    )

    @test length(data) == 100
end

@testset "Bulge: default mass & invalid bounds" begin
    @test Bulge(; NumSamples = 10).TotalMass == 8.57e9u"Msun"

    config = Bulge(; NumSamples = 10)
    @test_throws ErrorException generate(config; MaxRadius = 0.0u"kpc")
    @test_throws ErrorException generate(config; MaxHeight = 0.0u"kpc")
end

@testset "Bulge: collection propagation" begin
    config = Bulge(; NumSamples = 10, collection = GAS)
    particles = generate(config)
    @test all(==(GAS), particles.Collection)
end

@testset "Bulge: Base.show" begin
    config = Bulge(; NumSamples = 10)
    io = IOBuffer()
    show(io, config)
    @test !isempty(String(take!(io)))
end

using Random: seed!

@testset "Bulge: unit consistency" begin
    config = Bulge(;
        NumSamples  = 30,
        ScaleRadius = 0.075u"kpc",
        CutRadius   = 2.1u"kpc",
    )

    seed!(42)
    p_astro = generate(config)
    seed!(42)
    p_si = generate(config, uSI)

    # MaxRadius defaults to 5 * ScaleRadius = 0.375 kpc.
    for (a, b) in zip(p_astro, p_si)
        R_astro = sqrt(a.Pos.x^2 + a.Pos.y^2)
        R_si = sqrt(b.Pos.x^2 + b.Pos.y^2)
        @test ustrip(u"kpc", R_astro) ≈ ustrip(u"kpc", R_si)
        @test ustrip(u"kpc", a.Pos.z) ≈ ustrip(u"kpc", b.Pos.z)
        @test ustrip(u"kpc", R_astro) <= 0.375 + 1e-9
        @test abs(ustrip(u"kpc", a.Pos.z)) <= 0.375 + 1e-9
    end

    # Mixed units in the constructor are normalized to ScaleRadius's unit.
    config_pc = Bulge(;
        NumSamples  = 30,
        ScaleRadius = 75.0u"pc",
        CutRadius   = 2100.0u"pc",
    )
    p_pc = generate(config_pc)
    @test all(q -> sqrt(q.Pos.x^2 + q.Pos.y^2) <= 375.0u"pc" + 1e-6u"pc", p_pc)
end

@testset "Bulge: RotationCurve" begin
    config = Bulge(; NumSamples = 30)
    rc_r = [0.5, 1.0, 2.0, 4.0, 6.0] .* u"kpc"
    rc_v = [50.0, 100.0, 150.0, 200.0, 220.0] .* u"km/s"

    particles = generate(config; RotationCurve = (rc_r, rc_v), MaxRadius = 6.0u"kpc")
    @test length(particles) == 30
    @test all(p -> sqrt(p.Vel.x^2 + p.Vel.y^2) > 0u"km/s", particles)
end

@testset "Bulge: VelocityDispersion" begin
    seed!(7)
    config = Bulge(; NumSamples = 300)
    particles = generate(config; VelocityDispersion = 100.0u"km/s")

    vx = [ustrip(u"km/s", p.Vel.x) for p in particles]
    @test any(!iszero, vx)
    @test std(vx) ≈ 100.0 rtol = 0.2
    @test ustrip(u"km/s", norm(mean(particles.Vel))) < 20.0

    # Radial profile form
    rd = r -> 100.0u"km/s" * (r / 1.0u"kpc")
    p2 = generate(config; VelocityDispersion = rd)
    @test all(!iszero(norm(p.Vel)) for p in p2)
end

@testset "Bulge: MaxRadius warning" begin
    # MaxRadius < ScaleRadius → warn.
    config = Bulge(; NumSamples = 50, ScaleRadius = 1.0u"kpc")
    @test_logs (:warn, r"MaxRadius.*smaller than.*ScaleRadius") match_mode = :any generate(config; MaxRadius = 0.1u"kpc")
end
