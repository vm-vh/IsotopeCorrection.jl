var documenterSearchIndex = {"docs":
[{"location":"examples/","page":"Examples","title":"Examples","text":"EditURL = \"examples.jl\"","category":"page"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/#Example-1:-Basic-Usage","page":"Examples","title":"Example 1: Basic Usage","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"In the first example we will correct a measured response vector of the M-57 alanine fragment. First the chemical formula of the fragment and the response vector have to be defined, after which the uncorrected mass isotopomer distribution (MID) can be calculated.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"note: Labeled Carbon\nThe number after \"LabC\" in the chemical formula represents the number of possibly labeled carbon atoms in the fragment while the number behind \"C\" is the total number of carbon atoms, both from the alanine and from any derivative groups.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using IsotopeCorrection\n\nala_m57 = \"C11H26N1O2Si2LabC3\"\nresponse_vec = [2500000, 2000000, 2500000, 800000]\nMID = response_vec ./ sum(response_vec)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The function isotope_correction()` can be directly used on a single measurement vector. Here we are assuming that the tracer used for this experiment had an 80% purity.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"corr_resp, corr_MID, mean_enrich, residuum = isotope_correction(response_vec,\n                                                                ala_m57;\n                                                                tracer_purity = 0.8)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"corr_resp` is the corrected measurement vector,","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"corr_resp","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"corr_MID is the MID of the corrected measurement vector,","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"corr_MID","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"mean_enrich is the mean enrichment of the corrected measurement vector,","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"mean_enrich","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"and residuum is the residals of the optimization.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"residuum","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"The optimization can also be turned off, using the optimization = false keyword argument, though this isn't recommended as the optimization prevents negative values in the corrected response.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"To visualize the effects of isotope_correction(), the MID can be plotted before and after correction.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using ColorSchemes, CairoMakie\n\ncolors = colorschemes[:bamako][1:64:end];\ngroup_color = [PolyElement(color = color, strokecolor = :transparent)\n    for color in colors]\nlabels = [\"m+0\", \"m+1\", \"m+2\", \"m+3\"]\nMplus = [0, 1, 2, 3]\n\nfig = Figure();\nax = Axis(fig[1,1],\n    xticks = (1:2, [\"uncorrected\", \"corrected\"]),\n    title = \"Alanine (M-57) MID\")\nbarplot!(ax, repeat([1, 2], inner = length(corr_MID)),\n    cat(MID, corr_MID, dims = 1),\n    stack = repeat(Mplus, 2),\n    colormap = colors,\n    color = repeat(Mplus .+ 1, 2))\n\nLegend(fig[1,2], group_color, labels, framevisible = false)\nfig","category":"page"},{"location":"examples/#Example-2:-Usage-with-a-DataFrame","page":"Examples","title":"Example 2: Usage with a DataFrame","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"For this example we will make a very simple DataFrame with 2 samples for each of which 2 methionine fragments were measured, specifically M-159 and M-57. We need the chemical formulas of both fragments, as shown in the code below, as well as measured response vectors. For this example the measured response for sample 1 is all possible m+ being recorded with equal frequency while the measured response for sample 2 is simply a random vector.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"note: Labeled Carbon\nThe number after \"LabC\" in the chemical formula represents the number of possibly labeled carbon atoms in the fragment while the number behind \"C\" is the total number of carbon atoms, both from the methionine and from any derivative groups.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using IsotopeCorrection\nusing DataFrames, DataFramesMeta\n\ndf = DataFrame(\n    Sample = repeat([1, 2], inner = 11),\n    AminoAcidFragment = repeat(cat(repeat([\"Met_(M-159)\"], 5), repeat([\"Met_(M-57)\"], 6), dims = 1), 2),\n    Molecule = repeat(cat(repeat([\"C10H24N1S1Si1LabC4\"], 5), repeat([\"C13H30N1O2S1Si2LabC5\"], 6), dims = 1), 2),\n    MZ0 = repeat(cat(repeat([218], 5), repeat([320], 6), dims = 1), 2),\n    Mplus = repeat(cat(repeat(0:4, 1), repeat(0:5, 1), dims = 1), 2),\n    Response = cat(repeat([3e6/5], 5), repeat([3e6/6], 6), rand(0:3e6, 11), dims = 1))\n\n# Of course normally the first step would be to read in (and prepare) real data,\n# using CSV.jl or XLSX.jl, for example.\n\n# using CSV\n# df = CSV.read(\"DATA.csv\", DataFrame)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Using DataFramesMeta the DataFrame can be grouped based on the sample ID and amino acid fragment and each group corrected using the main function of this package: isotope_correction()","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Since this function outputs the corrected response, MID, mean enrichment and residuum, we index into element 1 to recieve only the corrected response. Additionally we are assuming that the tracer used for this experiment had a 99% purity.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"gdfs = groupby(df, [:Sample, :AminoAcidFragment])\ndf_corrected = @transform(gdfs,\n    :CorrectedResponse = isotope_correction(Vector(:Response), :Molecule[1]; tracer_purity = 0.99)[1])","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"Now, again using DataFramesMeta, both the corrected and uncorrected MIDs can be calculated.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"gdfs_corrected = groupby(df_corrected, [:Sample, :AminoAcidFragment])\ndf_corrected_MIDs = @transform(gdfs_corrected,\n    :MID = :Response ./ sum(:Response),\n    :CorrectedMID = :CorrectedResponse ./ sum(:CorrectedResponse))\nselect(df_corrected_MIDs, :Sample, :AminoAcidFragment, :Response, :CorrectedResponse, :MID, :CorrectedMID)\n# And of course the finished DataFrame can then be saved.\n# CSV.write(\"Corrected_MIDs.csv\", df_corrected_MIDs)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"To visualize what isotope_correction() does we can plot the MIDs before and after correction.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"using ColorSchemes, CairoMakie\n\ngdfs = @groupby(df_corrected_MIDs, :Sample, :AminoAcidFragment)\n\ncolors = colorschemes[:bamako][1:43:256];\ngroup_color = [PolyElement(color = color, strokecolor = :transparent)\n    for color in colors]\nlabels = [\"m+0\", \"m+1\", \"m+2\", \"m+3\", \"m+4\", \"m+5\"]\n\nfig = Figure();\n\nfunction myplot(group, fig, colors)\n    name = group.AminoAcidFragment[1]\n    sample = group.Sample[1]\n    ax = Axis(fig,\n        xticks = (1:2, [\"uncorrected\", \"corrected\"]),\n        title = \"Sample $sample $name MID\")\n    barplot!(ax, repeat([1, 2], inner = length(group.MID)),\n        cat(group.MID, group.CorrectedMID, dims = 1),\n        stack = repeat(group.Mplus, 2),\n        colormap = colors,\n        color = repeat(group.Mplus .+ 1, 2))\nend\nmyplot(gdfs[1], fig[1,1], colors[1:5])\nmyplot(gdfs[2], fig[1,2], colors)\nmyplot(gdfs[3], fig[2,1], colors[1:5])\nmyplot(gdfs[4], fig[2,2], colors)\n\nLegend(fig[1:2,3], group_color, labels, framevisible = false)\n\nfig\n# save(\"Methionine_MIDs.png\", fig)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"This page was generated using Literate.jl.","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"CurrentModule = IsotopeCorrection","category":"page"},{"location":"background/#Background","page":"Background","title":"Background","text":"","category":"section"},{"location":"background/","page":"Background","title":"Background","text":"Documentation for IsotopeCorrection.","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"Metabolic flux analysis (MFA) is a powerful method used to assess intracellular metabolic activity by estimating  the fluxes of different metabolic pathways. This approach is often more informative than simply measuring  metabolite concentrations as it helps to understand how cells respond to different conditions.  The stable isotope of carbon, carbon-13 (13C), is commonly used as a tracer in MFA to investigate metabolic  fluxes, though other isotopes like nitrogen-15 (15N) can be used. These isotopes do not significantly alter  the chemical properties of the molecules they label, making them ideal for tracking metabolic processes.  However, naturally occurring isotopes can affect the accuracy of the isotope tracing,  as their presence interferes with the signals obtained from the isotopically labeled substrates.","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"When conducting isotope labeling experiments, it is critical to distinguish between isotopes introduced into  the system experimentally and those naturally present at the start of the experiment. The natural abundance (NA)  of stable isotopes can falsify the results of isotope labeling studies if not corrected appropriately.{1} For example, carbon naturally exists as  approximately 98.93% 12C and 1.07% 13C.  Ignoring the NA of 13C could lead to errors in the interpretation of metabolic fluxes, especially  when using gas chromatography mass spectrometry (GCMS) as the derivatization necessary for gas chromatography  leads to more incorporated atoms.","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"While early correction methods were based on the assumption that mass isotopomer distributions (MIDs)  of labeled standards were simply shifted versions of the unlabeled standards, it is much more accurate  to account for the non-linear distribution of stable isotopes, or \"skew,\" caused by isotopic enrichment from  experimental tracers.  This package implements this “skewed” correction method, using multinomial probability theory to account  for the non-linear distribution of naturally occurring isotopes. A least-squares optimization function ensures  that fractional abundances remain positive after correction is used and corrects for noise,  similarly to IsoCor.{2} In addition to correcting for the natural abundances of all isotopes, correction for  the purity of the tracer substrate is also implemented.","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"{1}\tMidani, F. S.; Wynn, M. L.; Schnell, S. The Importance of Accurately Correcting for the Natural Abundance  of Stable Isotopes. Anal. Biochem. 2017, 520, 27–43. https://doi.org/10.1016/j.ab.2016.12.011.","category":"page"},{"location":"background/","page":"Background","title":"Background","text":"{2}\tMillard, P.; Letisse, F.; Sokol, S.; Portais, J.-C. IsoCor: Correcting MS Data in Isotope Labeling  Experiments. Bioinformatics 2012, 28 (9), 1294–1296. https://doi.org/10.1093/bioinformatics/bts127.","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"CurrentModule = IsotopeCorrection","category":"page"},{"location":"#IsotopeCorrection","page":"Reference","title":"IsotopeCorrection","text":"","category":"section"},{"location":"","page":"Reference","title":"Reference","text":"Documentation for IsotopeCorrection.","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"","category":"page"},{"location":"","page":"Reference","title":"Reference","text":"Modules = [IsotopeCorrection]","category":"page"},{"location":"#IsotopeCorrection.IsotopeCorrection","page":"Reference","title":"IsotopeCorrection.IsotopeCorrection","text":"A module for natural isotope abundance correction of fragments generated by GC-/LCMS in a 13C Isotope labeling experiment.\n\n\n\n\n\n","category":"module"},{"location":"#IsotopeCorrection.CM_conv-Tuple{String, String, Dict, Float64}","page":"Reference","title":"IsotopeCorrection.CM_conv","text":"CM_conv(\n    formula::String,\n    label::String,\n    isotopes::Dict,\n    tracer_purity::Float64\n) -> Matrix{Float64}\n\n\nOutputs the correction matrix constructed by iterative convolution of the correction vector for each column.\n\nArguments\n\nformula: full chemical formula of the fragment to be corrected e.g. \"C3H3O3LabC3\" (or \"C3H3O3\" in conjunction with label)\nlabel: labeled element and amount e.g. \"C3\", can be ignored if the formula contains LabC\nisotopes: dictionary containing the natural isotope abundances of elements, default is isotope_NAs\ntracer_purity::Float64: purity of the tracer, default is 1.0 (corresponding to 100%)\n\n\n\n\n\n","category":"method"},{"location":"#IsotopeCorrection.calculate_MDV-Tuple{String, String, Dict}","page":"Reference","title":"IsotopeCorrection.calculate_MDV","text":"calculate_MDV(\n    formula::String,\n    label::String,\n    isotopes::Dict\n) -> Tuple{Vector{Float64}, Int64}\n\n\nOutputs a tuple containing two elements:\n\nThe mass distribution vector (MDV) of a fragment, calculated by convolving the vectors of    natural isotope abundances n times for each element with n atoms.\nThe number of possible isotopologues (which is equal to the number of possibly labeled carbon atoms + 1).\n\nArguments\n\nformula: full chemical formula of the fragment to be corrected e.g. \"C3H3O3LabC3\" (or \"C3H3O3\" in conjunction with label)\nlabel: labeled element and amount e.g. \"C3\", can be ignored if the formula contains LabC\nisotopes: dictionary containing the natural isotope abundances of elements, default is isotope_NAs\n\n\n\n\n\n","category":"method"},{"location":"#IsotopeCorrection.fragment-Tuple{String, String}","page":"Reference","title":"IsotopeCorrection.fragment","text":"fragment(\n    formula::String,\n    label::String\n) -> Tuple{Dict{SubString{String}, Int64}, Int64}\n\n\nOutputs a tuple containing two elements:\n\nA dictionary of all elements in the fragment and their respective amounts.\nThe number of possible isotopologues (which is equal to the number of possibly labeled carbon atoms + 1).\n\nArguments\n\nformula::String: full chemical formula of the fragment to be corrected e.g. \"C3H3O3LabC3\" (or \"C3H3O3\" in conjunction with label)\nlabel::String: labeled element and amount e.g. \"C3\", can be ignored if the formula contains LabC\n\n\n\n\n\n","category":"method"},{"location":"#IsotopeCorrection.isotope_correction-Tuple{Vector, String}","page":"Reference","title":"IsotopeCorrection.isotope_correction","text":"isotope_correction(\n    response::Vector,\n    formula::String;\n    label,\n    isotopes,\n    tracer_purity,\n    optimization\n) -> Union{Tuple{Any, Any, Any}, NTuple{4, Any}}\n\n\nOutputs the corrected response measurements, mass isotopomer distribution (MID) and mean enrichment of a fragment. Additionally outputs the residuum when using L-BFGS-B for optimization. Only to be used with low resolution GC-/LCMS data.\n\nArguments\n\nresponse::Vector: the response / measurement vector obtained from LC-/GCMS, must have    the length n+1 if n is equal to the amount of possibly labeled atoms in the fragment\nformula::String: chemical formula of the fragment to be corrected.\nlabel::String = \"\": labeled element and amount e.g. \"C3\", can be ignored if the formula contains LabC\nisotopes::Dict = isotope_NAs: dictionary containing the natural isotope abundances of elements\ntracer_purity::Number = 1.0: purity of the tracer, the default represents 100%.\noptimization::Bool = true: whether to use the L-BFGS-B optimization algorithm.\n\nExamples\n\njulia> response_vec = [3500000, 1000000, 3500000, 800000]\njulia> pyruvate = \"C3H3O3\"\njulia> pyruvate_label = \"C3\"\njulia> corr_response, corr_MID, mean_enrichment, residuum = isotope_correction(response_vec, \n                                                                               pyruvate; \n                                                                               label = pyruvate_label)\n([3.6336363843460013e6, 844166.0884382844, 3.5694537799896887e6, 780365.311750165], \n[0.4116212229745574, 0.09562780668246623, 0.4043505664463884, 0.08840040389658813], \n0.5423825378162519, \n[-1.0054050521417098e-15, -2.645802768793973e-16, -8.466568860140713e-16, -1.852061938155781e-16])\n\njulia> pyruvate = \"C3H3O3LabC3\"\njulia> corr_response, corr_MID, mean_enrichment, residuum = isotope_correction(response_vec, \n                                                                               pyruvate; \n                                                                               tracer_purity = 0.99)\n([3.6336363843460013e6, 844166.0884382844, 3.5694537799896887e6, 780365.311750165], \n[0.4116212229745574, 0.09562780668246623, 0.4043505664463884, 0.08840040389658813], \n0.5423825378162519, \n[-1.0054050521417098e-15, -2.645802768793973e-16, -8.466568860140713e-16, -1.852061938155781e-16])\n\n\n\n\n\n","category":"method"}]
}
