@testset "Solar System" begin
    data = solarsystem(now())
    @test length(data) == 9
end