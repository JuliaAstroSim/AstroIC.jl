@info "Loading"
include("../../src/AstroIC.jl")
# include("AstroIC.jl\\src\\AstroIC.jl")

using Main.AstroIC
config = PlummerStarCluster()

@info "Generating"
data = generate(config)

using AstroPlot
plotrotationcurve(data)