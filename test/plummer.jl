@testset "Plummer" begin
    config = PlummerStarCluster(
        collection = STAR,
        NumSamples = 1000,
        VirialRadius = 0.010u"kpc",
        TotalMass = 1.0e5u"Msun",
        model = AstroIC.Newton(),
    )

    data = generate(
        config,
        MaxRadius = 0.050u"kpc",
    )

    @test length(data) == 1000
end