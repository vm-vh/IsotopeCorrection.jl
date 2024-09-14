using NaturalIsotopeCorrection
using Test

@testset "NaturalIsotopeCorrection.jl" begin
    include("../docs/src/examples.jl")
    #=
    pyruvate = "C3H3O3LabC3"
    response_vec = [3500000, 1000000, 3500000, 800000]

    corr_response, MID, mean_enrichment, residuum = isotope_correction(response_vec, pyruvate; tracer_purity = 0.99)
    @test corr_response ≈ [3.6690841286892043e6, 880069.4781954776, 3.548554810172248e6, 765218.8519553166]
    @test MID ≈ [0.41398107163956405, 0.09929783371601102, 0.4003818041674765, 0.08633929047694855]
    @test mean_enrichment ≈ 0.5397698283704524
    #@test residuum ≈ [-1.0858374563130465e-13, -2.566428685730154e-15, 1.1059455573558807e-13, -1.6554787924343889e-13]
    @test isapprox(residuum, [-1.0858374563130465e-13, -2.566428685730154e-15, 1.1059455573558807e-13, -1.6554787924343889e-13]; atol = 1e16)

    corr_response, MID, mean_enrichment = isotope_correction(response_vec, pyruvate; tracer_purity = 0.99, optimization = false)
    @test corr_response ≈ [3.6690841286882004e6, 880069.4781954774, 3.548554810173271e6, 765218.8519538215]
    @test MID ≈ [0.41398107163951964, 0.09929783371602752, 0.40038180416765856, 0.08633929047679423]
    @test mean_enrichment ≈ 0.5397698283704317
    =#
end
