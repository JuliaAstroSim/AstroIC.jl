# Tests for the extended Solar-System API in `src/data/Planets.jl` and
# `src/planets_extended.jl`. Network-dependent testsets are gated behind the
# `AstroIC_NET_TEST` environment variable so CI on offline / firewalled
# runners still passes.
#
#   - Always run: PLANET_TABLE integrity, lookup helpers, NASA consistency
#     checks for masses / radii / semi-major axes against published values.
#   - Run if `AstroIC_NET_TEST == "1"`: `horizons_state_vector`, `moon`,
#     `solar_system_full` (requires live access to ssd-api.jpl.nasa.gov).

# ---------------------------------------------------------------------------
# PLANET_TABLE integrity (pure, no network)
# ---------------------------------------------------------------------------

@testset "PLANET_TABLE: schema & lookup helpers" begin
    # Every entry must have a name, horizons_id, parent, and numeric fields.
    for p in PLANET_TABLE
        @test p.name isa AbstractString
        @test p.horizons_id isa AbstractString
        @test p.parent isa AbstractString
        @test p.a_AU isa Real
        @test p.e isa Real
        @test p.i_deg isa Real
        @test p.P_year isa Real
        @test p.r_eq_km isa Real
        @test p.M_kg isa Real
        @test p.albedo isa Real
        @test p.color isa Tuple{Float64, Float64, Float64}
        @test all(0.0 <= c <= 1.0 for c in p.color)
        @test 0.0 <= p.e < 1.0
        @test 0.0 <= p.albedo <= 1.0
        @test p.r_eq_km > 0
        @test p.M_kg > 0
    end

    # Names are unique.
    names = [p.name for p in PLANET_TABLE]
    @test length(names) == length(unique(names))

    # planet_names() must return exactly the table names.
    @test planet_names() == names
    @test length(planet_names()) == length(PLANET_TABLE)

    # Lookup round-trip.
    for p in PLANET_TABLE
        @test planet_info(p.name).name == p.name
        @test planet_horizons_id(p.name) == p.horizons_id
        @test planet_color(p.name) == p.color
    end

    # Unknown name throws.
    @test_throws ArgumentError planet_info("Pluton")
    @test_throws ArgumentError planet_horizons_id("Ceres")
    @test_throws ArgumentError planet_color("Sedna")
end

@testset "PLANET_TABLE: parent DAG & horizons-id sanity" begin
    # Sun is the root of the heliocentric family.
    @test planet_info("Sun").parent == ""

    # Every non-Sun body's parent must exist in the table.
    names = Set(planet_names())
    for p in PLANET_TABLE
        if p.parent != ""
            @test p.parent in names
        end
    end

    # No cycles: the ancestor chain from each body terminates at the Sun.
    function _ancestors(name)
        chain = String[]
        cur = name
        seen = Set{String}()
        while true
            par = planet_info(cur).parent
            par == "" && return chain
            @test par ∉ seen   # would mean a cycle
            push!(seen, par)
            push!(chain, par)
            cur = par
        end
    end
    # Every non-Sun body ultimately traces back to the Sun.
    for p in PLANET_TABLE
        if p.name == "Sun"
            continue
        end
        chain = _ancestors(p.name)
        @test "Sun" in chain
    end

    # Every horizons_id is a non-empty digit string (Horizons COMMAND ids
    # for major bodies use 3 digits like "199", "501", but Sun uses "10" and
    # Pluto uses "999"). This guards against typos like "x199" or " ".
    for p in PLANET_TABLE
        @test occursin(r"^\d+$", p.horizons_id)
        @test !isempty(p.horizons_id)
    end

    # Expected horizons_id mapping (spot-check). Cross-checked against
    # https://ssd.jpl.nasa.gov/horizons/.
    expected = Dict(
        "Sun"     => "10",
        "Mercury" => "199",
        "Venus"   => "299",
        "Earth"   => "399",
        "Moon"    => "301",
        "Mars"    => "499",
        "Jupiter" => "599",
        "Io"      => "501",
        "Europa"  => "502",
        "Ganymede"=> "503",
        "Callisto"=> "504",
        "Saturn"  => "699",
        "Titan"   => "606",
        "Uranus"  => "799",
        "Neptune" => "899",
        "Pluto"   => "999",
    )
    for (name, id) in expected
        @test planet_horizons_id(name) == id
    end
end

# ---------------------------------------------------------------------------
# NASA / published-value consistency (pure, no network)
# ---------------------------------------------------------------------------

# Tolerance helpers — PLANET_TABLE values are taken from NASA Planetary Fact
# Sheet & IAU 2015 B3. Tolerances are generous (1% for masses, 0.5% for
# semi-major axes, 2% for radii) so this test only catches real bugs
# (typos, swapped rows, factor-of-10 mistakes), not normal rounding.
@testset "PLANET_TABLE: NASA / IAU consistency" begin
    # Earth mass = 5.972e24 kg (CODATA / IAU 2015).
    @test isapprox(planet_info("Earth").M_kg,   5.972e24,  rtol = 0.01)
    # Jupiter mass = 1.898e27 kg.
    @test isapprox(planet_info("Jupiter").M_kg, 1.898e27,  rtol = 0.01)
    # Sun mass = 1.989e30 kg.
    @test isapprox(planet_info("Sun").M_kg,     1.989e30,  rtol = 0.01)
    # Moon mass = 7.342e22 kg.
    @test isapprox(planet_info("Moon").M_kg,    7.342e22,  rtol = 0.01)
    # Mercury semi-major axis = 0.38710 AU.
    @test isapprox(planet_info("Mercury").a_AU, 0.38710,  rtol = 0.005)
    # Neptune semi-major axis = 30.06896 AU.
    @test isapprox(planet_info("Neptune").a_AU, 30.06896, rtol = 0.005)
    # Pluto semi-major axis = 39.48169 AU.
    @test isapprox(planet_info("Pluto").a_AU,   39.48169, rtol = 0.005)
    # Earth equatorial radius = 6378.1 km.
    @test isapprox(planet_info("Earth").r_eq_km,   6378.1, rtol = 0.02)
    # Jupiter equatorial radius = 69911 km.
    @test isapprox(planet_info("Jupiter").r_eq_km, 69911.0, rtol = 0.02)
    # Saturn equatorial radius = 58232 km.
    @test isapprox(planet_info("Saturn").r_eq_km,  58232.0, rtol = 0.02)

    # Sun must be by far the heaviest body (at least 1000× Earth mass).
    m_sun = planet_info("Sun").M_kg
    @test m_sun > 1000 * planet_info("Earth").M_kg

    # Jupiter > Saturn > Neptune > Uranus > Earth > Venus > Mars > Mercury
    # (ordering check; uses table values directly — guards against swapped rows).
    order = ["Jupiter", "Saturn", "Neptune", "Uranus", "Earth",
             "Venus",   "Mars",   "Mercury"]
    for i in 1:(length(order) - 1)
        @test planet_info(order[i]).M_kg > planet_info(order[i + 1]).M_kg
    end
end

# ---------------------------------------------------------------------------
# Network-dependent tests: gated behind ASTROIC_NET_TEST=1
# ---------------------------------------------------------------------------

const _NET_TEST = get(ENV, "ASTROIC_NET_TEST", "0") == "1"

if _NET_TEST
    @testset "horizons_state_vector (network)" begin
        # Sun heliocentric position at J2000 is at the SSB (within ~R_sun).
        sv = horizons_state_vector("10", DateTime(2000, 1, 1, 12, 0);
                                   center = "500@10")
        @test length(sv) == 6
        x, y, z, vx, vy, vz = sv
        # Position in km: |r| of the Sun w.r.t. the SSB is ≲ a few × 10^5 km.
        @test sqrt(x^2 + y^2 + z^2) < 1.0e6
        # Velocity in km/s: barycentric velocity of the Sun is ≲ 0.015 km/s.
        @test sqrt(vx^2 + vy^2 + vz^2) < 0.05

        # Earth heliocentric distance at J2000 should be ~1 AU = 1.496e8 km.
        sv = horizons_state_vector("399", DateTime(2000, 1, 1, 12, 0);
                                   center = "500@10")
        x, y, z, _, _, _ = sv
        r_au = sqrt(x^2 + y^2 + z^2) / 1.496e8
        @test 0.98 < r_au < 1.02

        # Cache: identical call must hit the cache (returns identical tuple).
        sv_a = horizons_state_vector("10", DateTime(2025, 1, 1))
        sv_b = horizons_state_vector("10", DateTime(2025, 1, 1))
        @test sv_a === sv_b
    end

    @testset "moon (network)" begin
        # Default Moon relative to Earth.
        m = moon(DateTime(2025, 1, 1, 0, 0))
        @test m isa Star
        @test m.Mass > 0
        @test m.Pos isa PVector
        @test m.Vel isa PVector
        # Geocentric distance ~ 3.84e5 km (lunar semi-major axis).
        r_km = ustrip(u"km", norm(m.Pos))
        @test 3.5e5 < r_km < 4.2e5
        # Mass should match planet_info("Moon").M_kg.
        @test isapprox(ustrip(u"kg", m.Mass), planet_info("Moon").M_kg,
                       rtol = 1e-12)
    end

    @testset "solar_system_full (network)" begin
        # Default: 9 (sun + 8 planets) + Pluto + 6 moons = 16.
        bodies = solar_system_full(DateTime(2025, 1, 1))
        @test length(bodies) == 16
        # All particles should have non-zero mass and finite Pos.
        for b in bodies
            @test b isa Star
            @test b.Mass > 0
            @test all(isfinite, (ustrip(u"m", b.Pos.x),
                                 ustrip(u"m", b.Pos.y),
                                 ustrip(u"m", b.Pos.z)))
        end

        # Toggle flags.
        b_no_pluto = solar_system_full(DateTime(2025, 1, 1);
                                       include_pluto = false)
        @test length(b_no_pluto) == 15
        b_no_moons = solar_system_full(DateTime(2025, 1, 1);
                                       include_moons = false)
        @test length(b_no_moons) == 10   # sun + 8 planets + Pluto
        b_minimal = solar_system_full(DateTime(2025, 1, 1);
                                      include_pluto = false,
                                      include_moons = false)
        @test length(b_minimal) == 9
    end

    @testset "moon: unknown parent (network)" begin
        # Unknown parent raises before any network call.
        @test_throws ErrorException moon(DateTime(2025, 1, 1);
                                        body = "Moon", parent = "Pluto")
    end
else
    @testset "Network tests skipped (set ASTROIC_NET_TEST=1 to enable)" begin
        @test true
    end
end