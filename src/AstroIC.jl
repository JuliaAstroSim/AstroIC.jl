module AstroIC

__precompile__(true)

using Unitful, UnitfulAstro
using Distributions
using Random
using BangBang
using Dates

using PhysicalParticles
using AstroIO

import Base: show
import Unitful: Units
import PhysicalConstants: CODATA2018

using AstroLib
import AstroLib: planets

export
    show,

    # Physics
    vmean,

    GravModel,
        MOND,
        Newton,

    PlummerStarCluster,
    GasCloud,

    solarsystem,

    generate

abstract type InitialConditionConfig end

include("Traits.jl")
include("physics.jl")

include("plummer.jl")
include("gascloud.jl")
include("solarsystem.jl")

end # module
