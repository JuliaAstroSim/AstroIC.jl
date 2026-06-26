@testset "Bulge" begin
    config = Bulge(;
        collection = STAR,
        NumSamples = 100,
        TotalMass = 8.57u"Msun",
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

@testset "Bulge: Base.show" begin
    config = Bulge(; NumSamples = 10)
    io = IOBuffer()
    show(io, config)
    @test !isempty(String(take!(io)))
end

@testset "Bulge: MaxRadius warning" begin
    # MaxRadius < ScaleRadius → warn.
    config = Bulge(; NumSamples = 50, ScaleRadius = 1.0u"kpc")
    @test_logs (:warn, r"MaxRadius.*smaller than.*ScaleRadius") match_mode = :any generate(config; MaxRadius = 0.1u"kpc")
end
