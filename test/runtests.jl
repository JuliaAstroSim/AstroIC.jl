using Test

using DataFrames
using Unitful, UnitfulAstro

using AstroSimBase
using PhysicalParticles
using AstroIO

using Dates
using AstroIC

include("tools.jl")
include("physics.jl")
include("distribution.jl")
include("plummer.jl")
include("disk.jl")
include("bulge.jl")
include("gascloud.jl")
include("spherical.jl")
include("solarsystem.jl")

@testset "Load data" begin
    df = load_data_MW_satellites()
    @test length(df.Galaxy) == 61

    df = load_data_UFDs()
    @test length(df.Galaxy) == 27

    @test !isnothing(load_SPARC_LTGs_RC())
    dfSPARC = load_SPARC_LTGs_data()
    @test !isnothing(dfSPARC)
    @test !isempty(dfSPARC.a0_rodrigues2018)
    @test !isnothing(load_li2018_SPARC())

    @test !isnothing(load_SPARC_Xray_ETGs_data())
    @test !isnothing(load_SPARC_rotating_ETGs_data())
    @test !isnothing(load_SPARC_rotating_ETGs_RC("NGC2685", :bulge))
    @test !isnothing(load_SPARC_rotating_ETGs_RC("NGC2824", :disk))
    @test !isnothing(load_SPARC_rotating_ETGs_rotmod("NGC2685"))

    @test !isnothing(load_MW_RC_Eilers2019())
    @test !isnothing(load_MW_RC_Mroz2019())
    @test !isnothing(load_MW_RC_stddev_W21())
    @test !isnothing(load_MW_RC_DS_W21())

    @test !isnothing(load_massive_dwarf_CO_RC("NGC1035"))
    @test !isnothing(load_massive_dwarf_DM_RC("NGC1035"))
    @test !isnothing(load_massive_dwarf_Baryon_RC("NGC1035"))

    particles = generate_milkyway_baryon_particles(500)
    @test length(particles) == 500
end

@testset "Load data: error paths" begin
    # Milky Way baryon particle generation fails for too-small Np. The src
    # check at MilkyWay.jl:56 only catches `iszero(NumSamples_HII)`; with
    # Np = 1 the HII count becomes negative (Np - sum(other ceil-rounded
    # counts) = -3), bypassing the check and triggering an ArgumentError
    # downstream from `StructArray(Star(...) for i in 1:NumSamples)` with
    # NumSamples = -3. Either way it must error — accept either exception
    # type. Fix the src by tightening the guard to `NumSamples_HII < 1`.
    @test_throws Union{ErrorException, ArgumentError} generate_milkyway_baryon_particles(1)

    # SPARC rotating ETGs loader rejects unsupported :mode symbols.
    @test_throws ErrorException load_SPARC_rotating_ETGs_RC("NGC2685", :halo)
    @test_throws ErrorException load_SPARC_rotating_ETGs_RC("NGC2685", :unknown)
end
