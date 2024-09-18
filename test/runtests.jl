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
    @test isapprox(residuum, [-1.0858374563130465e-13, -2.566428685730154e-15, 1.1059455573558807e-13, -1.6554787924343889e-13]; atol = 1e16)

    corr_response, MID, mean_enrichment = isotope_correction(response_vec, pyruvate; tracer_purity = 0.99, optimization = false)
    @test corr_response ≈ [3.6690841286882004e6, 880069.4781954774, 3.548554810173271e6, 765218.8519538215]
    @test MID ≈ [0.41398107163951964, 0.09929783371602752, 0.40038180416765856, 0.08633929047679423]
    @test mean_enrichment ≈ 0.5397698283704317
    
    ala_m57 = "C11H26N1O2Si2LabC3"
    response_vec = [2500000, 2000000, 2500000, 800000]
    corr_response, corr_MID, mean_enrichment, residuum = isotope_correction(response_vec, ala_m57; tracer_purity = 0.8)

    @test corr_response ≈ [3.7494356377578774e6, 1.582246275076472e6, 3.110758584864337e6, 392021.5450371306]
    @test corr_MID ≈ [0.42441018135800024, 0.17909933478943207, 0.3521163563572243, 0.04437412749534335]
    @test mean_enrichment ≈ 0.5041136074974777
    @test isapprox(residuum, [7.161631798132872e-13, 8.285188904175392e-13, -1.0657673462843284e-12, 6.363440591555375e-13]; atol = 1e16)

    corr_response, corr_MID, mean_enrichment = isotope_correction(response_vec, pyruvate; tracer_purity = 0.8, optimization = false)
    @test corr_response ≈ [2.788815753935856e6, 1.7402842780064365e6, 2.8703208164970432e6, 961526.0409580052]
    @test corr_MID ≈ [0.333552621590313, 0.20814440051201863, 0.3433009268527882, 0.11500205104488019]
    @test mean_enrichment ≈ 0.5599381018380588
    =#
end
