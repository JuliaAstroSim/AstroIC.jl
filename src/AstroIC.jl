module AstroIC

__precompile__(true)

using Reexport
using Unitful, UnitfulAstro
using Distributions
using Random
using BangBang
using Dates
using StructArrays

@reexport using PhysicalParticles
@reexport using AstroIO

import Base: show
import Unitful: Units
import PhysicalConstants: CODATA2018
import PhysicalParticles: rotate, rotate_x, rotate_y, rotate_z

using AstroLib
import AstroLib: planets

export
    # Tools
    setpos,
    setvel,
    addpos,
    addvel,

    # Physics
    vmean,

    GravModel,
        MOND,
        Newton,

    PlummerStarCluster,
    GasCloud,

    solarsystem,

    helio2xyz,
    gridpoints,

    generate

abstract type InitialConditionConfig end

_common_keywords = """
## Common keywords
- `constants`
"""

include("Traits.jl")
include("Tools.jl")
include("physics.jl")

include("plummer.jl")
include("gascloud.jl")
include("solarsystem.jl")


"""
    function generate(::InitialConditionConfig, units; kw...)

Generate initial conditions in `units`

# Usable IC configs
- `PlummerStarCluster`
- `GasCloud`
- `ExponentialDisk`

$_common_keywords
"""
generate

end # module
