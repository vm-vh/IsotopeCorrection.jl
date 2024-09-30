using IsotopeCorrection
using Test, DataFrames, DataFramesMeta

@testset "IsotopeCorrection.jl" begin
    # include("../docs/src/examples.jl")

    ala_m57 = "C11H26N1O2Si2LabC3"
    response_vec = [2500000, 2000000, 2500000, 800000]

    corr_response, corr_MID, mean_enrichment, residuum = isotope_correction(response_vec, ala_m57; tracer_purity = 0.8)
    @test corr_response ≈ [3.026602632345966e6, 819243.5134320166, 3.6669247890182203e6, 587930.8175325772]
    @test corr_MID ≈ [0.37362227679544946, 0.10113241278096696, 0.4526675467300172, 0.07257776369356636]
    @test mean_enrichment ≈ 0.5560501993304251
    @test isapprox(residuum, [-2.022641591536693e-13, -2.999336291582156e-13, 3.7879754717533407e-13, 2.285023816885092e-13]; atol = 1e16)

    corr_response, corr_MID, mean_enrichment = isotope_correction(response_vec, ala_m57; tracer_purity = 0.8, optimization = false)
    @test corr_response ≈ [3.026602632344788e6, 819243.5134261829, 3.666924789023547e6, 587930.8175357856]
    @test corr_MID ≈ [0.3736222767952338, 0.1011324127802278, 0.4526675467305896, 0.07257776369394879]
    @test mean_enrichment ≈ 0.5560501993308133

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

    @test df_corrected_MIDs.CorrectedResponse[1:11] ≈ [253750.6191521022, 207468.52394624625, 191847.0947432641, 
        195070.04435999502, 206287.42079718813, 238333.70861583538, 174867.32075263848, 
        155844.14228048082, 162681.3506664206, 163463.03446118528, 173492.2420107591
    ]
    @test df_corrected_MIDs.CorrectedMID[1:11] ≈ [0.24065337153407296, 0.1967601101494616, 0.1819449754379082, 
        0.18500157366077136, 0.19563996921778581, 0.22301653203627414, 0.16362898755370225, 
        0.14582838638902995, 0.1522261826214522, 0.1529576293398783, 0.162342282059663
    ]
end
