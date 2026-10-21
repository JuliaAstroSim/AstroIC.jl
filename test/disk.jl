@testset "Disk" begin
    config = ExponentialDisc(
        collection = STAR,
        NumSamples = 100,
        TotalMass = 1.0e8u"Msun",
        ScaleRadius = 2.0u"kpc",
        ScaleHeight = 0.02u"kpc",
    )

    data = generate(
        config,
        MaxRadius = 4.0u"kpc",
    )

    @test length(data) == 100
end

@testset "Disk: invalid sampling bounds" begin
    config = ExponentialDisc(; NumSamples = 10)
    @test_throws ErrorException generate(config; MaxRadius = 0.0u"kpc")
    @test_throws ErrorException generate(config; MaxHeight = 0.0u"kpc")
end

@testset "Disk: collection propagation" begin
    config = ExponentialDisc(; NumSamples = 10, collection = GAS)
    particles = generate(config)
    @test all(==(GAS), particles.Collection)
end

@testset "Disk: Base.show" begin
    config = ExponentialDisc(; NumSamples = 10)
    io = IOBuffer()
    show(io, config)
    @test !isempty(String(take!(io)))
end

@testset "Disk: MaxRadius / MaxHeight warnings" begin
    config = ExponentialDisc(;
        NumSamples   = 50,
        ScaleRadius  = 5.0u"kpc",
        ScaleHeight  = 0.5u"kpc",
    )

    # MaxRadius < ScaleRadius → warn.
    @test_logs (:warn, r"MaxRadius.*smaller than.*ScaleRadius") match_mode = :any generate(config; MaxRadius = 0.5u"kpc")

    # MaxHeight < ScaleHeight → warn.
    @test_logs (:warn, r"MaxHeight.*smaller than.*ScaleHeight") match_mode = :any generate(config; MaxRadius = 50.0u"kpc", MaxHeight = 0.05u"kpc")
end

using Random: seed!

@testset "Disk: unit consistency" begin
    config = ExponentialDisc(;
        NumSamples   = 30,
        ScaleRadius  = 2.0u"kpc",
        ScaleHeight  = 0.02u"kpc",
    )

    seed!(42)
    p_astro = generate(config; MaxRadius = 6.0u"kpc")
    seed!(42)
    p_si = generate(config, uSI; MaxRadius = 6.0u"kpc")

    # The same config must describe the same physical system in either unit
    # system, and the sampled radii must respect MaxRadius/MaxHeight.
    for (a, b) in zip(p_astro, p_si)
        R_astro = sqrt(a.Pos.x^2 + a.Pos.y^2)
        R_si = sqrt(b.Pos.x^2 + b.Pos.y^2)
        @test ustrip(u"kpc", R_astro) ≈ ustrip(u"kpc", R_si)
        @test ustrip(u"kpc", a.Pos.z) ≈ ustrip(u"kpc", b.Pos.z)
        @test ustrip(u"kpc", R_astro) <= 6.0 + 1e-9
        @test abs(ustrip(u"kpc", a.Pos.z)) <= 6.0 + 1e-9
    end

    # Length arguments may be given in a different unit than the output unit.
    config_pc = ExponentialDisc(;
        NumSamples   = 30,
        ScaleRadius  = 2000.0u"pc",
        ScaleHeight  = 20.0u"pc",
    )
    p_pc = generate(config_pc)
    @test all(
        q -> sqrt(q.Pos.x^2 + q.Pos.y^2) <= 10000.0u"pc" + 1e-6u"pc",
        p_pc,
    )
    @test all(q -> abs(q.Pos.z) <= 10000.0u"pc" + 1e-6u"pc", p_pc)
end

@testset "Disk: with RotationCurve" begin
    # Providing a rotation curve switches on the Spline1D + rotational_velocity path.
    config = ExponentialDisc(;
        NumSamples   = 50,
        ScaleRadius  = 2.0u"kpc",
        ScaleHeight  = 0.02u"kpc",
    )
    rc_r = [0.5, 1.0, 2.0, 4.0, 6.0] .* u"kpc"
    rc_v = [50.0, 100.0, 150.0, 200.0, 220.0] .* u"km/s"

    particles = generate(config; RotationCurve = (rc_r, rc_v), MaxRadius = 6.0u"kpc")
    @test length(particles) == 50

    # Velocities should be non-zero in the x/y plane (rotational).
    for p in particles
        v_mag = sqrt(p.Vel.x^2 + p.Vel.y^2)
        @test v_mag > 0u"km/s"
    end
end
