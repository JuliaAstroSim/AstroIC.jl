@testset "GasCloud: gridpoints" begin
    R = 1.0u"kpc"
    Nx, Ny, Nz = 3, 3, 3
    x, y, z = gridpoints(R, Nx, Ny, Nz)
    @test size(x) == (Nx, Ny, Nz)
    @test size(y) == (Nx, Ny, Nz)
    @test size(z) == (Nx, Ny, Nz)

    # x varies along dim 1, y along dim 2, z along dim 3.
    @test all(x[:, 1, 1] .≈ collect(-R:(2*R/(Nx-1)):R))
    @test all(y[1, :, 1] .≈ collect(-R:(2*R/(Ny-1)):R))
    @test all(z[1, 1, :] .≈ collect(-R:(2*R/(Nz-1)):R))

    # Non-cubic grid.
    x2, _, _ = gridpoints(R, 5, 3, 2)
    @test size(x2) == (5, 3, 2)
end

@testset "GasCloud: generate (known buggy — see src/gascloud.jl:130)" begin
    # `generate(::GasCloud, ...)` references bare `SPHGas` at gascloud.jl:130,
    # but `SPHGas` is not exported by PhysicalParticles, so
    # `@reexport using PhysicalParticles` does not bring it into the AstroIC
    # namespace and the call raises `UndefVarError`.
    #
    # Fix in src: either `import PhysicalParticles: SPHGas` (or `using
    # PhysicalParticles: SPHGas`) at the top of gascloud.jl, or qualify the
    # constructor call as `PhysicalParticles.SPHGas(units; id)`.
    #
    # This test documents the bug; it should be flipped to a real assertion
    # (length > 0, all inside sphere, type <: SPHGas) once src is fixed.
    config = GasCloud(;
        Radius       = 1.0u"kpc",
        rho0         = 1.0u"Msun/kpc^3",
        T            = 100u"K",
        ParticleMass = Constant(uSI).m_p,
        Nx           = 5,
        Ny           = 5,
        Nz           = 5,
    )
    @test_throws UndefVarError generate(config, uSI)
end