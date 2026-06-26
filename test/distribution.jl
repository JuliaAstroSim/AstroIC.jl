@testset "Distribution: exponential family" begin
    # These helpers are not exported from AstroIC, so we access them via
    # the module.
    beta = 2.0

    # exponential_pdf: x > 0 branch
    @test AstroIC.exponential_pdf(beta, 1.0) ≈ exp(-0.5) / 2.0
    # x = 0 and x < 0 → else branch
    @test AstroIC.exponential_pdf(beta,  0.0) == 0.0
    @test AstroIC.exponential_pdf(beta, -1.0) == 0.0

    # exponential_cdf: x > 0 branch
    @test AstroIC.exponential_cdf(beta, 1.0) ≈ 1 - exp(-0.5)
    @test AstroIC.exponential_cdf(beta,  0.0) == 0.0
    @test AstroIC.exponential_cdf(beta, -1.0) == 0.0

    # exponential_cdf_inv: closed-form identity
    u = 0.7
    @test AstroIC.exponential_cdf_inv(beta, u) ≈ -beta * log(1 - u)

    # exponential_pdf_hole: x > 0 and else branches
    Rm = 0.5
    @test AstroIC.exponential_pdf_hole(beta, Rm,  1.0) ≈ exp(-1.0/beta - Rm/1.0) / beta
    @test AstroIC.exponential_pdf_hole(beta, Rm,  0.0) == 0.0
    @test AstroIC.exponential_pdf_hole(beta, Rm, -1.0) == 0.0
end

@testset "Distribution: sech2 family" begin
    # sech2_pdf: identity  sech(x)^2 = 1 / cosh(x)^2
    @test AstroIC.sech2_pdf(0.0) ≈ 1.0
    @test AstroIC.sech2_pdf(1.0) ≈ 1 / cosh(1.0)^2

    # sech2_cdf: cdf(0) = 0.5 (tanh(0) = 0)
    @test AstroIC.sech2_cdf(0.0)  ≈ 0.5
    @test AstroIC.sech2_cdf(Inf)  ≈ 1.0
    @test AstroIC.sech2_cdf(-Inf) ≈ 0.0

    # sech2_cdf_inv: round-trip with sech2_cdf
    u = 0.3
    @test AstroIC.sech2_cdf_inv(u) ≈ atanh(2*u - 1)
    @test AstroIC.sech2_cdf(AstroIC.sech2_cdf_inv(0.5)) ≈ 0.5
end

@testset "Distribution: rejection_sampling" begin
    # Array-based dispatch: x, pdf_array, rMax, NumSamples.
    x = collect(0.0:0.01:5.0)
    pdf_array = exp.(-x ./ 1.0)   # exponential pdf, scale 1
    rMax = 5.0
    samples = AstroIC.rejection_sampling(x, pdf_array, rMax, 50)
    @test length(samples) == 50
    @test all(s -> 0 <= s <= rMax, samples)

    # Function-based dispatch.
    f = t -> exp(-t)
    samples2 = AstroIC.rejection_sampling(f, 5.0, 50)
    @test length(samples2) == 50
    @test all(s -> 0 <= s <= 5.0, samples2)
end
