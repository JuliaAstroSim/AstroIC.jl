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
