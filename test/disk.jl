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