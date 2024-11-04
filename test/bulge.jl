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