module AstroIC

__precompile__(true)

using PrecompileTools
using DocStringExtensions
using Reexport
using Unitful, UnitfulAstro
using Distributions
using Random
using BangBang
using Dates
using StructArrays
using Dierckx

@reexport using PhysicalParticles
@reexport using AstroIO

import Base: show
import Unitful: Units
import PhysicalConstants: CODATA2018
import PhysicalParticles: rotate, rotate_x, rotate_y, rotate_z
using AstroSimBase

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

    PlummerStarCluster,
    GasCloud,
    ExponentialDisk,
    Bulge,

    solarsystem,

    helio2xyz,
    gridpoints,

    generate

abstract type InitialConditionConfig end

_common_keywords = """
## Common keywords
- `constants`
"""

include("Tools.jl")
include("physics.jl")

include("plummer.jl")
include("disk.jl")
include("bulge.jl")
include("gascloud.jl")
include("solarsystem.jl")

include("precompile.jl")

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
