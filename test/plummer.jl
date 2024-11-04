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