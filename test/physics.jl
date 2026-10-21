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

@testset "Physics: vmean2" begin
    v2 = vmean2(300u"K", 1.0u"u")
    @test v2 ≈ vmean(300u"K", 1.0u"u")^2
    @test unit(v2) == getuVel(uAstro)^2

    v2_si = vmean2(300u"K", 1.0u"u", uSI)
    @test v2_si ≈ vmean(300u"K", 1.0u"u", uSI)^2
    @test unit(v2_si) == getuVel(uSI)^2
end
