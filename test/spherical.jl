# using AstroPlot

@testset "SphericalSystem" begin
    r_s = 2.0u"kpc"
    MaxRadius = 3*r_s
    NumSamples = 100

    r = collect(LinRange(0.001u"kpc", MaxRadius, 100))
    dr = r[2] - r[1]

    shell_mass_func = x->1.0e8u"Msun/kpc^3"/(x/r_s)/(1+x/r_s)^2 * 4π * x^2

    config = SphericalSystem(STAR, NumSamples, shell_mass_func)
    particles = generate(config; MaxRadius)
    # plot_makie(particles)


    shell_mass = shell_mass_func.(r)

    config = SphericalSystem(STAR, NumSamples, Dict("r"=>r, "m"=>shell_mass))
    particles = generate(config; MaxRadius)
    # plot_makie(particles)
end

@testset "SphericalSystem: Base.show" begin
    config = SphericalSystem(STAR, 100, Dict("r"=>[1.0u"kpc"], "m"=>[1.0e8u"Msun/kpc"]))
    io = IOBuffer()
    show(io, config)
    @test !isempty(String(take!(io)))
end

@testset "SphericalSystem: error paths" begin
    r_s = 2.0u"kpc"
    MaxRadius = 3 * r_s
    NumSamples = 50
    shell_mass_func = x -> 1.0e8u"Msun/kpc^3" / (x/r_s) / (1 + x/r_s)^2 * 4π * x^2

    # MaxRadius = 0 → error
    @test_throws ErrorException generate(SphericalSystem(STAR, NumSamples, shell_mass_func); MaxRadius = 0.0u"kpc")

    # Dict with missing keys → error
    bad_dict = Dict("a" => 1, "b" => 2)
    @test_throws ErrorException generate(SphericalSystem(STAR, NumSamples, bad_dict); MaxRadius = MaxRadius)

    # DataFrame with missing properties → error
    bad_df = DataFrame(x = [1, 2, 3])
    @test_throws ErrorException generate(SphericalSystem(STAR, NumSamples, bad_df); MaxRadius = MaxRadius)

    # Unsupported type → error
    @test_throws ErrorException generate(SphericalSystem(STAR, NumSamples, "not a valid type"); MaxRadius = MaxRadius)
end

@testset "SphericalSystem: DataFrame mass_shell (known buggy — see src/spherical.jl:65)" begin
    # The DataFrame `mass_shell` branch in `generate` has a bug at
    # spherical.jl:65: the TotalMass computation uses Dict-style indexing
    # `config.mass_shell["r"]` / `["m"]`, but DataFrame does not support
    # `df[column]` with a String — only `df[!, column]` or `df.column`.
    # The branching at line 50 (Dict vs DataFrame) uses `.r` / `.m` correctly,
    # but the TotalMass calculation at line 65-67 re-uses Dict syntax in the
    # shared `Dict || DataFrame` branch.
    #
    # Fix in src: split the TotalMass branch into Dict and DataFrame, or use
    # `getproperty(config.mass_shell, :r)` / `getproperty(..., :m)`.
    r_s = 2.0u"kpc"
    MaxRadius = 3 * r_s
    NumSamples = 50
    r = collect(LinRange(0.001u"kpc", MaxRadius, 50))
    shell_mass = 1.0e8u"Msun/kpc^3" ./ (r ./ r_s) ./ (1 .+ r ./ r_s).^2 .* 4π .* r.^2

    df = DataFrame(r = r, m = shell_mass)
    config = SphericalSystem(STAR, NumSamples, df)
    @test_throws ArgumentError generate(config; MaxRadius)
end
