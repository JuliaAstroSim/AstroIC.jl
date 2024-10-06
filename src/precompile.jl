@setup_workload begin
    @compile_workload begin
        NumSamples = 5
        data = generate(PlummerStarCluster(; NumSamples))
        data = generate(ExponentialDisc(; NumSamples))
        data = generate(ExponentialDisc(; NumSamples, HoleRadius = 0.2u"kpc"))
        data = generate(Bulge(; NumSamples))
        data = generate(Bulge(; NumSamples), RotationCurve = ([0.0, 5.0, 10.0]*u"kpc", [0.0, 200.0, 230.0]*u"km/s"))
        p = sum(data.Pos)
        v = sum(data.Vel)
        pos = PVector(100.0, 100.0, 100.0, u"kpc")
        vel = PVector(100.0, 100.0, 100.0, u"kpc/Gyr")
        setpos(data, pos)
        setvel(data, vel)

        solarsystem(now())
    end
end