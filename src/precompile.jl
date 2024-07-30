@setup_workload begin
    @compile_workload begin
        NumSamples = 5
        data = generate(PlummerStarCluster(; NumSamples))
        data = generate(ExponentialDisk(; NumSamples))
        p = sum(data.Pos)
        v = sum(data.Vel)
        pos = PVector(100.0, 100.0, 100.0, u"kpc")
        vel = PVector(100.0, 100.0, 100.0, u"kpc/Gyr")
        setpos(data, pos)
        setvel(data, vel)

        solarsystem(now())
    end
end