@testset "Physics: vmean" begin
    # vmean — Maxwell-Boltzmann thermal speed, default uAstro units.
    v = vmean(300u"K", 1.0u"u")
    @test v isa Quantity
    @test v > 0u"m/s"

    # vmean with explicit units argument (uSI).
    vSI = vmean(300u"K", 1.0u"u", uSI)
    @test vSI isa Quantity
    @test vSI > 0u"m/s"
end

@testset "Physics: vmean2 (known buggy — see src/physics.jl:6)" begin
    # `vmean2` is currently broken: the body computes
    #   8 * k_B * T / pi / m   (dimensions of velocity^2)
    # but then calls
    #   uconvert(getuVel(units), ...)   (expects velocity)
    # The conversion raises a DimensionError.
    #
    # The function also appears to be unused elsewhere in the codebase, so
    # this test is here as a regression marker until src is fixed (either
    # drop the `uconvert` and return the squared-speed Quantity as-is, or
    # remove the function entirely if it's not needed).
    @test_throws Unitful.DimensionError AstroIC.vmean2(300u"K", 1.0u"u")
end
