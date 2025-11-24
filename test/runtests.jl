using Test

using Unitful, UnitfulAstro

using AstroSimBase
using PhysicalParticles
using AstroIO

using Dates
using AstroIC

include("tools.jl")
include("plummer.jl")
include("disk.jl")
include("bulge.jl")
include("spherical.jl")
include("solarsystem.jl")

@testset "Load data" begin
    df = load_data_MW_satellites()
    @test length(df.Galaxy) == 61

    @test !isnothing(load_SPARC_LTGs_RC())
    @test !isnothing(load_SPARC_LTGs_data())
    @test !isnothing(load_li2018_SPARC())
    
    @test !isnothing(load_SPARC_Xray_ETGs_data())
    @test !isnothing(load_SPARC_rotating_ETGs_data())
    @test !isnothing(load_SPARC_rotating_ETGs_RC("NGC2685", :bulge))
    @test !isnothing(load_SPARC_rotating_ETGs_rotmod("NGC2685"))

    @test !isnothing(load_MW_RC_Eilers2019())
    @test !isnothing(load_MW_RC_Mroz2019())
    @test !isnothing(load_MW_RC_stddev_W21())
    @test !isnothing(load_MW_RC_DS_W21())
end