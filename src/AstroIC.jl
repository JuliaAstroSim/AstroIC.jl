module AstroIC

__precompile__(true)

using PrecompileTools
using DocStringExtensions
using Reexport
using Unitful, UnitfulAstro
using DataFrames
using Distributions
using Random
using BangBang
using Dates
using StructArrays
using Dierckx
using QuadGK
using Optim
using CSV
using UnitfulAngles
using Measurements

@reexport using PhysicalParticles
@reexport using AstroIO

import Base: show
import Unitful: Units
import Distributions: pdf
import PhysicalConstants: CODATA2018
import PhysicalParticles: rotate, rotate_x, rotate_y, rotate_z
@reexport using AstroSimBase

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
    ExponentialDisc,
    Bulge,
    SphericalSystem,

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
include("distribution.jl")

include("plummer.jl")
include("disk.jl")
include("bulge.jl")
include("spherical.jl")
include("gascloud.jl")
include("solarsystem.jl")

include("data/MW_satellites.jl")
export load_data_MW_satellites

include("data/UFDs.jl")
export load_data_UFDs

include("data/SPARC.jl")
export load_SPARC_LTGs_RC, load_SPARC_LTGs_data, load_li2018_SPARC
export load_SPARC_Xray_ETGs_data, load_SPARC_rotating_ETGs_data, load_SPARC_rotating_ETGs_RC, load_SPARC_rotating_ETGs_rotmod

include("data/MilkyWay.jl")
export load_MW_RC_Eilers2019, load_MW_RC_Mroz2019, load_MW_RC_stddev_W21, load_MW_RC_DS_W21

include("precompile.jl")

"""
    function generate(::InitialConditionConfig, units; kw...)

Generate initial conditions in `units`

# Usable IC configs
- `PlummerStarCluster`
- `GasCloud`
- `ExponentialDisc`

$_common_keywords
"""
generate

end # module
