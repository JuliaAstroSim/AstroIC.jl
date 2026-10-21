@testset "Plummer" begin
    config = PlummerStarCluster(
        collection = STAR,
        NumSamples = 100,
        VirialRadius = 0.010u"kpc",
        TotalMass = 1.0e5u"Msun",
        model = AstroSimBase.Newton(),
    )

    data = generate(
        config,
        MaxRadius = 0.050u"kpc",
    )

    @test length(data) == 100
end

@testset "Plummer: Base.show" begin
    config = PlummerStarCluster(; NumSamples = 10)
    io = IOBuffer()
    show(io, config)
    @test !isempty(String(take!(io)))
end

@testset "Plummer: MaxRadius = 0 disables truncation" begin
    config = PlummerStarCluster(; NumSamples = 10, VirialRadius = 0.01u"kpc")
    particles = generate(config; MaxRadius = 0.0u"kpc")
    @test length(particles) == 10
end

@testset "Plummer: collection propagation" begin
    config = PlummerStarCluster(; NumSamples = 10, collection = BLACKHOLE)
    particles = generate(config)
    @test all(==(BLACKHOLE), particles.Collection)
end

@testset "Plummer: MOND model" begin
    # Exercises the MOND1983Milgrom dispatch of plummer_vel_sigma2.
    config = PlummerStarCluster(;
        NumSamples    = 50,
        VirialRadius  = 0.01u"kpc",
        TotalMass     = 1.0e5u"Msun",
        model         = MOND1983Milgrom(),
    )
    particles = generate(config)
    @test length(particles) == 50
end

@testset "Plummer: unsupported gravity model" begin
    config = PlummerStarCluster(; NumSamples = 10, model = QUMOND())
    @test_throws ErrorException generate(config)
end

@testset "Plummer: MaxRadius warning" begin
    # MaxRadius < VirialRadius → warn.
    config = PlummerStarCluster(; NumSamples = 50, VirialRadius = 1.0u"kpc")
    @test_logs (:warn, r"MaxRadius.*smaller than.*VirialRadius") match_mode = :any generate(config; MaxRadius = 0.1u"kpc")
end
