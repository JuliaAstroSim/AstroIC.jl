var documenterSearchIndex = {"docs":
[{"location":"manual/plummer/#Plummer","page":"Plummer","title":"Plummer","text":"","category":"section"},{"location":"manual/plummer/","page":"Plummer","title":"Plummer","text":"PlummerStarCluster\ngenerate(::PlummerStarCluster)","category":"page"},{"location":"manual/plummer/#AstroIC.generate-Tuple{PlummerStarCluster}","page":"Plummer","title":"AstroIC.generate","text":"function generate(config::PlummerStarCluster, units = uAstro; kw...)\n\nGenerate initial conditions of Plummer model\n\nKeywords\n\nMaxRadius: resample particles outside the interested radius. Default is 5 * VirialRadius. Set to zero to avoid cutting off.\n\nCommon keywords\n\nconstants\n\n\n\n\n\n","category":"method"},{"location":"manual/plummer/","page":"Plummer","title":"Plummer","text":"using AstroIC\nusing UnitfulAstro\n\n# config\nconfig = PlummerStarCluster(\n    collection = STAR,\n    NumSamples = 1000,\n    VirialRadius = 0.010u\"kpc\",\n    TotalMass = 1.0e5u\"Msun\",\n    model = AstroIC.Newton(),\n)\n\n# generate\nparticles = generate(\n    config,\n    MaxRadius = 0.050u\"kpc\",\n)","category":"page"},{"location":"manual/plummer/","page":"Plummer","title":"Plummer","text":"tip: Use [AstroPlot.jl](https://github.com/JuliaAstroSim/AstroPlot.jl) to visualize\n","category":"page"},{"location":"lib/Methods/#Methods","page":"Methods","title":"Methods","text":"","category":"section"},{"location":"lib/Methods/#Index","page":"Methods","title":"Index","text":"","category":"section"},{"location":"lib/Methods/","page":"Methods","title":"Methods","text":"Pages = [\"Methods.md\"]","category":"page"},{"location":"lib/Methods/","page":"Methods","title":"Methods","text":"generate\naddpos\naddvel\nsetpos\nsetvel\nsolarsystem\nhelio2xyz","category":"page"},{"location":"lib/Methods/#AstroIC.generate","page":"Methods","title":"AstroIC.generate","text":"function generate(::InitialConditionConfig, units; kw...)\n\nGenerate initial conditions in units\n\nUsable IC configs\n\nPlummerStarCluster\nGasCloud\nExponentialDisk\n\nCommon keywords\n\nconstants\n\n\n\n\n\n","category":"function"},{"location":"lib/Methods/#AstroIC.addpos","page":"Methods","title":"AstroIC.addpos","text":"function addpos(data::Array, pos::AbstractPoint)\nfunction addpos(data::StructArray, pos::AbstractPoint)\n\nAdd pos to :Pos of all particles\n\n\n\n\n\n","category":"function"},{"location":"lib/Methods/#AstroIC.addvel","page":"Methods","title":"AstroIC.addvel","text":"function addvel(data::Array, vel::AbstractPoint)\nfunction addvel(data::StructArray, vel::AbstractPoint)\n\nAdd vel to :Vel of all particles\n\n\n\n\n\n","category":"function"},{"location":"lib/Methods/#AstroIC.setpos","page":"Methods","title":"AstroIC.setpos","text":"function setpos(data::Union{Array, StructArray}, pos::AbstractPoint)\n\nSet system center (middle value) to pos\n\n\n\n\n\n","category":"function"},{"location":"lib/Methods/#AstroIC.setvel","page":"Methods","title":"AstroIC.setvel","text":"function setvel(data::Union{Array, StructArray}, vel::AbstractPoint)\n\nSet system velocity (mass weighted average) to vel\n\n\n\n\n\n","category":"function"},{"location":"lib/Methods/#AstroIC.solarsystem","page":"Methods","title":"AstroIC.solarsystem","text":"function solarsystem(date::DateTime)\nfunction solarsystem(date::Real)\n\nGenerate initial conditions of Solar System at desinated date\n\n\n\n\n\n","category":"function"},{"location":"lib/Methods/#AstroIC.helio2xyz","page":"Methods","title":"AstroIC.helio2xyz","text":"function helio2xyz(jd, num)\n\nConvert heliocentric coordinates (in unit AU) to Cartesian coordinates in unit m\n\nReturns a PhysicalParticles::PVector\n\n\n\n\n\n","category":"function"},{"location":"lib/Types/#Types","page":"Types","title":"Types","text":"","category":"section"},{"location":"lib/Types/#Index","page":"Types","title":"Index","text":"","category":"section"},{"location":"lib/Types/","page":"Types","title":"Types","text":"Pages = [\"Types.md\"]","category":"page"},{"location":"lib/Types/","page":"Types","title":"Types","text":"PlummerStarCluster","category":"page"},{"location":"lib/Types/#AstroIC.PlummerStarCluster","page":"Types","title":"AstroIC.PlummerStarCluster","text":"struct PlummerStarCluster\n\nfields\n\ncollection particle type\nNumSamples amount of particles\nVirialRadius\nTotalMass mass are seperated equally to all particles\nG Newtonian constant of gravitation\nmodel gravity model [ Newton | MOND ]\n\n\n\n\n\n","category":"type"},{"location":"manual/guide/#Package-Guide","page":"Package Guide","title":"Package Guide","text":"","category":"section"},{"location":"manual/guide/#Installation","page":"Package Guide","title":"Installation","text":"","category":"section"},{"location":"manual/guide/","page":"Package Guide","title":"Package Guide","text":"From the Julia REPL, type ] to enter the Pkg REPL mode and run","category":"page"},{"location":"manual/guide/","page":"Package Guide","title":"Package Guide","text":"pkg> add AstroIC","category":"page"},{"location":"manual/guide/","page":"Package Guide","title":"Package Guide","text":"or add from git repository","category":"page"},{"location":"manual/guide/","page":"Package Guide","title":"Package Guide","text":"pkg> add https://github.com/JuliaAstroSim/AstroIC.jl","category":"page"},{"location":"manual/guide/","page":"Package Guide","title":"Package Guide","text":"Test the package by","category":"page"},{"location":"manual/guide/","page":"Package Guide","title":"Package Guide","text":"pkg> test AstroIC","category":"page"},{"location":"manual/guide/#Basic-Usage","page":"Package Guide","title":"Basic Usage","text":"","category":"section"},{"location":"manual/guide/","page":"Package Guide","title":"Package Guide","text":"using AstroIC\n\nconfig = PlummerStarCluster()\ndata = generate(config)","category":"page"},{"location":"#AstroIC.jl","page":"Home","title":"AstroIC.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Generates initial conditions for astrophysical simulations","category":"page"},{"location":"","page":"Home","title":"Home","text":"Source code: https://github.com/JuliaAstroSim/AstroIC.jl","category":"page"},{"location":"#Package-Feature","page":"Home","title":"Package Feature","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Plummer model\nSolar system","category":"page"},{"location":"#Manual-Outline","page":"Home","title":"Manual Outline","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"manual/solarsystem/#Solar-System","page":"Solar System","title":"Solar System","text":"","category":"section"},{"location":"manual/solarsystem/","page":"Solar System","title":"Solar System","text":"solarsystem","category":"page"},{"location":"manual/solarsystem/","page":"Solar System","title":"Solar System","text":"using AstroIC, Dates\n\nsolarsystem(now())","category":"page"},{"location":"manual/solarsystem/","page":"Solar System","title":"Solar System","text":"tip: Use [AstroPlot.jl](https://github.com/JuliaAstroSim/AstroPlot.jl) to visualize\n","category":"page"},{"location":"manual/tools/#Tools","page":"Tools","title":"Tools","text":"","category":"section"},{"location":"manual/tools/","page":"Tools","title":"Tools","text":"Here are some tools to manipulate initial conditions","category":"page"},{"location":"manual/tools/#Example:-collisions","page":"Tools","title":"Example: collisions","text":"","category":"section"},{"location":"manual/tools/","page":"Tools","title":"Tools","text":"config = PlummerStarCluster(NumSamples = 1000)\ngalaxy1 = generate(config);\ngalaxy2 = generate(config);\ngalaxy3 = generate(config);\ngalaxy4 = generate(config);\n\nsetpos(galaxy1, PVector(-0.1, 0.0, 0.0, u\"kpc\"))\nsetpos(galaxy2, PVector(+0.1, 0.0, 0.0, u\"kpc\"))\nsetpos(galaxy3, PVector(0.0, +0.1, 0.0, u\"kpc\"))\nsetpos(galaxy4, PVector(0.0, -0.1, 0.0, u\"kpc\"))\n\nsetvel(galaxy1, PVector(+0.4, +0.8, 0.0, u\"kpc/Gyr\"))\nsetvel(galaxy2, PVector(-0.4, -1.0, -0.2, u\"kpc/Gyr\"))\nsetvel(galaxy3, PVector(+1.0, -0.4, 0.4, u\"kpc/Gyr\"))\nsetvel(galaxy4, PVector(-1.0, +0.4, -0.2, u\"kpc/Gyr\"))\n\ndata = deepcopy(galaxy1);\nappend!(data, galaxy2);\nappend!(data, galaxy3);\nappend!(data, galaxy4);\ndata.ID .= 1:4000;","category":"page"}]
}
