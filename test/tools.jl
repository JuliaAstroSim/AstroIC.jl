@testset "Tools" begin
    NumSamples = 100
    data = generate(PlummerStarCluster(; NumSamples))
    p = sum(data.Pos)
    v = sum(data.Vel)
    pos = PVector(100.0, 100.0, 100.0, u"kpc")
    vel = PVector(100.0, 100.0, 100.0, u"kpc/Gyr")
    setpos(data, pos)
    setvel(data, vel)

    # estimate the error
    @test 1 - norm(mean(data.Pos)) / norm(pos) < 0.0001
    @test 1 - norm(mean(data.Vel)) / norm(vel) < 0.0001
end

@testset "Tools: addpos / addvel on Vector" begin
    # The StructArray branch is already covered by the "Tools" testset above
    # (the result of `generate` is a StructArray). This exercises the
    # `addpos(data::Array, pos)` and `addvel(data::Array, vel)` methods which
    # build a fresh Star via `setproperties!!` for every particle.

    n = 50
    arr = [Star(uAstro, id = i) for i in 1:n]
    @test length(arr) == n

    shift_pos = PVector(1.0, 2.0, 3.0, u"kpc")
    shift_vel = PVector(10.0, 20.0, 30.0, u"kpc/Gyr")

    addpos(arr, shift_pos)
    @test all(p -> p.Pos == shift_pos, arr)

    addvel(arr, shift_vel)
    @test all(p -> p.Vel == shift_vel, arr)

    # Round-trip: median should now match the shift.
    @test median(arr, :Pos) ≈ shift_pos
end

@testset "Tools: rotational_velocity family" begin
    # Place test point in the disk plane (z = 0) to make the cross product
    # easy to reason about. These helpers are not exported from AstroIC, so
    # we access them via the module.
    x, y = 1.0u"kpc", 0.0u"kpc"
    v    = 200.0u"km/s"

    # rot_ratio = 1.0 → pure rotational component (no random kick).
    v_pure = AstroIC.rotational_velocity(x, y, v, 1.0)
    @test v_pure isa PVector
    # At (1, 0) on a +z disk the tangential direction is ±y, no z-component.
    @test iszero(ustrip(u"km/s", v_pure.z))
    @test iszero(ustrip(u"km/s", v_pure.x))
    @test abs(ustrip(u"km/s", v_pure.y)) ≈ ustrip(u"km/s", v) atol = 1e-6

    # rot_ratio = 0.0 → fully random direction, magnitude still equals |v|.
    v_random = AstroIC.rotational_velocity(x, y, v, 0.0)
    @test v_random isa PVector
    @test norm(v_random) ≈ v atol = 1e-6 * u"km/s"

    # Default rot_ratio (no kwarg) is also 1.0.
    v_default = AstroIC.rotational_velocity(x, y, v)
    @test v_default isa PVector
    @test iszero(ustrip(u"km/s", v_default.z))
end

@testset "Tools: rotational_velocity_acc + freefall_velocity_acc" begin
    # `a` is acceleration magnitude; result magnitude should be sqrt(a * r).
    a = 0.01u"kpc/Gyr^2"
    x, y, z = 1.0u"kpc", 2.0u"kpc", 0.5u"kpc"

    vec_acc = AstroIC.rotational_velocity_acc(x, y, z, a, 1.0)
    @test vec_acc isa PVector
    r = sqrt(ustrip(u"kpc", x)^2 + ustrip(u"kpc", y)^2 + ustrip(u"kpc", z)^2)
    v_expected = sqrt(ustrip(u"kpc/Gyr^2", a) * r)
    @test norm(vec_acc) ≈ v_expected * u"kpc/Gyr" rtol = 1e-6

    # freefall_velocity_acc: pure radial infall, magnitude also sqrt(a*r).
    vec_free = AstroIC.freefall_velocity_acc(x, y, z, a)
    @test vec_free isa PVector
    @test norm(vec_free) ≈ v_expected * u"kpc/Gyr" rtol = 1e-6
end
