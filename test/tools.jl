@testset "Tools" begin
    NumSamples = 100
    data = generate(PlummerStarCluster(; NumSamples))
    p = sum(data.Pos)
    v = sum(data.Vel)
    pos = PVector(100.0, 100.0, 100.0, u"kpc")
    vel = PVector(100.0, 100.0, 100.0, u"kpc/Gyr")
    setpos(data, pos)
    setvel(data, vel)

    # estimate the error
    @test 1 - norm(mean(data.Pos)) / norm(pos) < 0.0001
    @test 1 - norm(mean(data.Vel)) / norm(vel) < 0.0001
end