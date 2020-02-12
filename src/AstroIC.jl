module AstroIC

__precompile__(true)

using Unitful, UnitfulAstro
using Distributions
using Random

using PhysicalParticles
using AstroIO

import Base: show
import Unitful: Units

export
    show,

    PlummerStarCluster,

    generate

abstract type InitialConditionConfig{T} end

include("Traits.jl")
include("plummer.jl")

end # module
