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
end