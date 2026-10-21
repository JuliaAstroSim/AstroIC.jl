using AstroLib: helio, jdcnv

@testset "Solar System" begin
    jd = jdcnv(DateTime(2024, 1, 1))
    data = solarsystem(jd)
    @test length(data) == 9

    @testset "heliocentric Cartesian conversion" begin
        # `helio` returns (radius, longitude, latitude); the Cartesian
        # components must be the standard spherical-to-Cartesian conversion.
        for i in 1:8
            r, longitude, latitude = helio(jd, i, true)
            p = data.Pos[i]
            @test ustrip(u"AU", p.x) ≈ r * cos(latitude) * cos(longitude)
            @test ustrip(u"AU", p.y) ≈ r * cos(latitude) * sin(longitude)
            @test ustrip(u"AU", p.z) ≈ r * sin(latitude)
            @test ustrip(u"AU", norm(p)) ≈ r
        end
    end

    @testset "Earth position & velocity" begin
        p = data.Pos[3]
        v = data.Vel[3]

        # Earth stays close to the ecliptic plane.
        @test abs(ustrip(u"AU", p.z)) / ustrip(u"AU", norm(p)) < 0.001

        # Nearly circular orbit: the radial velocity component is small.
        radial_fraction = abs(p * v) / (norm(p) * norm(v))
        @test radial_fraction < 0.05

        # Earth's orbital speed is ~29.8 km/s.
        @test 27.0u"km/s" < norm(v) < 33.0u"km/s"
    end
end
