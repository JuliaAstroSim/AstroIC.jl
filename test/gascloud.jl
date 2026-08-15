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

@testset "GasCloud: generate unit consistency" begin
    # Bug history: `generate(::GasCloud, uSI)` previously crashed with
    #   MethodError: Cannot `convert` an object of type Star to an object
    #   of type Star
    # because `pos` was in kpc (from `gridpoints(R)`) and `mass` was in Msun
    # (from `rho0 = Msun/kpc^3`), while `Star(uSI)` parameterizes Pos in m
    # and Mass in kg. `setproperties!!` then triggers a `convert` that
    # BangBang cannot reconcile, and you get the cryptic "convert Star to
    # Star" MethodError. The fix in src/gascloud.jl converts Pos to the
    # target length unit and Mass to the target mass unit before assignment.

    for units in (uSI, uAstro)
        config = GasCloud(;
            Radius       = 1.0u"kpc",
            rho0         = 1.0u"Msun/kpc^3",
            T            = 100u"K",
            ParticleMass = Constant(units).m_p,
            Nx           = 5,
            Ny           = 5,
            Nz           = 5,
        )
        data = generate(config, units)
        @test length(data) > 0

        # Every particle must live inside the sphere and carry finite Pos.
        u_len = getuLength(units)
        R_stripped = ustrip(u_len, 1.0u"kpc")   # = 1.0 (kpc) or 3.086e19 (m)
        for p in data
            r = ustrip(u_len, norm(p.Pos))
            @test r ≤ R_stripped
            @test all(isfinite, (ustrip(u_len, p.Pos.x),
                                 ustrip(u_len, p.Pos.y),
                                 ustrip(u_len, p.Pos.z)))
            # Mass must be in the Star's native mass unit (no leftover Msun).
            @test unit(p.Mass) == getuMass(units)
            @test ustrip(getuMass(units), p.Mass) > 0
        end
    end
end