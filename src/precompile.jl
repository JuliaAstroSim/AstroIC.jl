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

        load_data_MW_satellites()
        load_data_UFDs()

        load_SPARC_LTGs_RC()
        load_SPARC_LTGs_data()
        load_li2018_SPARC()
        
        load_SPARC_Xray_ETGs_data()
        load_SPARC_rotating_ETGs_data()
        load_SPARC_rotating_ETGs_RC("NGC2685", :bulge)
        load_SPARC_rotating_ETGs_rotmod("NGC2685")

        load_MW_RC_Eilers2019()
        load_MW_RC_Mroz2019()
        load_MW_RC_stddev_W21()
        load_MW_RC_DS_W21()

        generate_milkyway_baryon_particles(500)

        load_massive_dwarf_CO_RC("NGC1035")
        load_massive_dwarf_DM_RC("NGC1035")
    end
end