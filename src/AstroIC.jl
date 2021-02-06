module AstroIC

__precompile__(true)

using Unitful, UnitfulAstro
using Distributions
using Random
using BangBang

using PhysicalParticles
using AstroIO

import Base: show
import Unitful: Units
import PhysicalConstants: CODATA2018

export
    show,

    # Physics
    vmean,

    GravModel,
        MOND,
        Newton,

    PlummerStarCluster,
    GasCloud,

    generate

abstract type InitialConditionConfig end

include("Traits.jl")
include("physics.jl")

include("plummer.jl")
include("gascloud.jl")

end # module
