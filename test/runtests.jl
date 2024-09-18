using NaturalIsotopeCorrection
using Test, DataFrames, DataFramesMeta

@testset "NaturalIsotopeCorrection.jl" begin
    # include("../docs/src/examples.jl")

    ala_m57 = "C11H26N1O2Si2LabC3"
    response_vec = [2500000, 2000000, 2500000, 800000]

    corr_response, corr_MID, mean_enrichment, residuum = isotope_correction(response_vec, ala_m57; tracer_purity = 0.8)
    @test corr_response ≈ [3.7494356377578774e6, 1.582246275076472e6, 3.110758584864337e6, 392021.5450371306]
    @test corr_MID ≈ [0.42441018135800024, 0.17909933478943207, 0.3521163563572243, 0.04437412749534335]
    @test mean_enrichment ≈ 0.5041136074974777
    @test isapprox(residuum, [7.161631798132872e-13, 8.285188904175392e-13, -1.0657673462843284e-12, 6.363440591555375e-13]; atol = 1e16)

    corr_response, corr_MID, mean_enrichment = isotope_correction(response_vec, ala_m57; tracer_purity = 0.8, optimization = false)
    @test corr_response ≈ [3.7494356377637773e6, 1.5822462750899016e6, 3.1107585848453073e6, 392021.54504765704]
    @test corr_MID ≈ [0.424410181358148, 0.17909933479073276, 0.35211635635463884, 0.044374127496480495]
    @test mean_enrichment ≈ 0.504113607497363

    df = DataFrame(
        Sample = repeat([1, 2], inner = 11), 
        AminoAcidFragment = repeat(cat(repeat(["Met_(M-159)"], 5), repeat(["Met_(M-57)"], 6), dims = 1), 2),
        Molecule = repeat(cat(repeat(["C10H24N1S1Si1LabC4"], 5), repeat(["C13H30N1O2S1Si2LabC5"], 6), dims = 1), 2),
        MZ0 = repeat(cat(repeat([218], 5), repeat([320], 6), dims = 1), 2),
        Mplus = repeat(cat(repeat(0:4, 1), repeat(0:5, 1), dims = 1), 2),
        Response = cat(repeat([1e6/5], 5), repeat([1e6/6], 6), rand(0:1e6, 11), dims = 1)
    )

    gdfs = groupby(df, [:Sample, :AminoAcidFragment])
    df_corrected = @transform(gdfs, 
        :CorrectedResponse = isotope_correction(Vector(:Response), :Molecule[1]; tracer_purity = 0.99)[1]
    )

    gdfs_corrected = groupby(df_corrected, [:Sample, :AminoAcidFragment])
    df_corrected_MIDs = @transform(gdfs_corrected, 
        :MID = :Response ./ sum(:Response),
        :CorrectedMID = :CorrectedResponse ./ sum(:CorrectedResponse)
    )

    @test df_corrected_MIDs.CorrectedResponse[1:11] ≈ [256350.5953947725, 209316.21271342717, 
        213614.23098928804, 213323.38585836615, 215543.78052654164, 
        240727.618300901, 176403.71542171348, 168185.8260919769, 
        173128.23316293195, 173266.76220269877, 174721.82172768374
    ]
    @test df_corrected_MIDs.CorrectedMID[1:11] ≈ [0.2313324103459418, 0.18888828378538802, 
        0.19276684285771886, 0.19250438235876846, 0.19450808065218303, 
        0.21757070311023, 0.15943447065382058, 0.1520071053511907, 
        0.15647407507022204, 0.1565992782388322, 0.15791436757570465
    ]
end
