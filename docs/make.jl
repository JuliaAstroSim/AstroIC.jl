"""
Compile with:
julia --project=docs/ --color=yes docs/make.jl

Generate key:
DocumenterTools.genkeys(user="JuliaAstroSim", repo="git@github.com:JuliaAstroSim/AstroIC.jl.git")
"""

using Documenter

using AstroIC

# The DOCSARGS environment variable can be used to pass additional arguments to make.jl.
# This is useful on CI, if you need to change the behavior of the build slightly but you
# can not change the .travis.yml or make.jl scripts any more (e.g. for a tag build).
if haskey(ENV, "DOCSARGS")
    for arg in split(ENV["DOCSARGS"])
        (arg in ARGS) || push!(ARGS, arg)
    end
end

makedocs(
    modules = [AstroIC],
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        prettyurls = !("local" in ARGS),
        canonical = "https://juliaastrosim.github.io/AstroIC.jl/dev/",
        assets = ["assets/alpha_small.ico"],
        analytics = "UA-153693590-1",
        highlights = ["llvm", "yaml"],
    ),
    clean = false,
    sitename = "AstroIC.jl",
    authors = "islent",
    linkcheck = !("skiplinks" in ARGS),
    linkcheck_ignore = [
        # Repository no longer under the JuliaAstroSim GitHub org
        "https://github.com/JuliaAstroSim/ISLENT",
        # IOPscience / Oxford Academic URLs that return bot-protection 302 / 403
        # but are still the canonical references for the cited papers.
        "https://iopscience.iop.org/article/10.3847/1538-4357/aaf648",  # Eilers+ 2019, ApJ 871, 120
        "https://iopscience.iop.org/article/10.3847/2041-8213/aaf73f",  # Mróz+ 2019, ApJ 870, L10
        # GitHub URLs in src/index.md (Package ecosystem list) — curl times out
        # intermittently from CI but the repos do exist.
        "https://github.com/JuliaAstroSim/PhysicalParticles.jl",
        "https://github.com/JuliaAstroSim/AstroIO.jl",
        "https://github.com/JuliaAstroSim/AstroIC.jl",
        "https://github.com/JuliaAstroSim/ParallelOperations.jl",
        "https://github.com/JuliaAstroSim/PhysicalTrees.jl",
        "https://github.com/JuliaAstroSim/PhysicalMeshes.jl",
        "https://github.com/JuliaAstroSim/AstroPlot.jl",
        "https://github.com/JuliaAstroSim/AstroNbodySim.jl",
        "https://github.com/JuliaAstroSim/WaveDM.jl",
        "https://github.com/JuliaAstro/AstroLib.jl",
    ],
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "manual/guide.md",
            "manual/plummer.md",
            "manual/disk.md",
            "manual/bulge.md",
            "manual/spherical.md",
            "manual/gascloud.md",
            "manual/solarsystem.md",
            "manual/tools.md",
            "manual/data.md",
        ],
        "Library" => Any[
            "lib/Types.md",
            "lib/Methods.md",
        ],
        #"contributing.md",
    ],
    #strict = !("strict=false" in ARGS),
    #doctest = ("doctest=only" in ARGS) ? :only : true,
)

deploydocs(
    repo = "github.com/JuliaAstroSim/AstroIC.jl.git",
    target = "build",
)