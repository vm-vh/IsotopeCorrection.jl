# # Examples

# ## Example 1: Basic Usage
#=
In the first example we will correct a measured response vector of the M-57 alanine fragment.
First the chemical formula of the fragment and the response vector have to be defined, after 
which the uncorrected mass isotopomer distribution (MID) can be calculated.
=#

#md # !!! note "Labeled Carbon"
#md #     The number after "LabC" in the chemical formula represents the number of possibly
#md #     labeled carbon atoms in the fragment while the number behind "C" is the total number of 
#md #     carbon atoms, both from the alanine and from any derivative groups. 

using IsotopeCorrection

ala_m57 = "C11H26N1O2Si2LabC3"
response_vec = [2500000, 2000000, 2500000, 800000]
MID = response_vec ./ sum(response_vec)

# The function `isotope_correction()`` can be directly used on a single measurement vector.
# Here we are assuming that the tracer used for this experiment had an 80% purity.
corr_resp, corr_MID, mean_enrich, residuum = isotope_correction(response_vec, 
                                                                ala_m57; 
                                                                tracer_purity = 0.8)

# `corr_resp`` is the corrected measurement vector,
corr_resp

# `corr_MID` is the MID of the corrected measurement vector,
corr_MID

# `mean_enrich` is the mean enrichment of the corrected measurement vector,
mean_enrich

# and `residuum` is the residals of the optimization.
residuum

# The optimization can also be turned off, using the `optimization = false` keyword argument,
# though this isn't recommended as the optimization prevents negative values in the corrected response.


# To visualize the effects of `isotope_correction()`, 
# the MID can be plotted before and after correction.
using ColorSchemes, CairoMakie

colors = colorschemes[:bamako][1:64:end];
group_color = [PolyElement(color = color, strokecolor = :transparent)
    for color in colors]
labels = ["m+0", "m+1", "m+2", "m+3"]
Mplus = [0, 1, 2, 3]

fig = Figure();
ax = Axis(fig[1,1], 
    xticks = (1:2, ["uncorrected", "corrected"]),
    title = "Alanine (M-57) MID")
barplot!(ax, repeat([1, 2], inner = length(corr_MID)), 
    cat(MID, corr_MID, dims = 1),
    stack = repeat(Mplus, 2),
    colormap = colors,
    color = repeat(Mplus .+ 1, 2))

Legend(fig[1,2], group_color, labels, framevisible = false)
fig



# ## Example 2: Usage with a DataFrame

#=
For this example we will make a very simple DataFrame with 2 samples for each of which
2 methionine fragments were measured, specifically M-159 and M-57.
We need the chemical formulas of both fragments, as shown in the code below, as well as
measured response vectors. For this example the measured response for sample 1 is 
all possible m+ being recorded with equal frequency while the measured response for 
sample 2 is simply a random vector.
=#

#md # !!! note "Labeled Carbon"
#md #     The number after "LabC" in the chemical formula represents the number of possibly
#md #     labeled carbon atoms in the fragment while the number behind "C" is the total number of 
#md #     carbon atoms, both from the methionine and from any derivative groups. 

using IsotopeCorrection
using DataFrames, DataFramesMeta

df = DataFrame(
    Sample = repeat([1, 2], inner = 11), 
    AminoAcidFragment = repeat(cat(repeat(["Met_(M-159)"], 5), repeat(["Met_(M-57)"], 6), dims = 1), 2),
    Molecule = repeat(cat(repeat(["C10H24N1S1Si1LabC4"], 5), repeat(["C13H30N1O2S1Si2LabC5"], 6), dims = 1), 2),
    MZ0 = repeat(cat(repeat([218], 5), repeat([320], 6), dims = 1), 2),
    Mplus = repeat(cat(repeat(0:4, 1), repeat(0:5, 1), dims = 1), 2),
    Response = cat(repeat([1e6/5], 5), repeat([1e6/6], 6), rand(0:1e6, 11), dims = 1))

## Of course normally the first step would be to read in (and prepare) real data,
## using CSV.jl or XLSX.jl, for example.

## using CSV
## df = CSV.read("DATA.csv", DataFrame)

# Using DataFramesMeta the DataFrame can be grouped based on the sample ID and amino acid 
# fragment and each group corrected using the main function of this package: 
# `isotope_correction()`

# Since this function outputs the corrected response, MID, mean enrichment and residuum,
# we index into element 1 to recieve only the corrected response.
# Additionally we are assuming that the tracer used for this experiment had a 99% purity.
gdfs = groupby(df, [:Sample, :AminoAcidFragment])
df_corrected = @transform(gdfs, 
    :CorrectedResponse = isotope_correction(Vector(:Response), :Molecule[1]; tracer_purity = 0.99)[1])

# Now, again using DataFramesMeta, both the corrected and uncorrected MIDs can be calculated.
gdfs_corrected = groupby(df_corrected, [:Sample, :AminoAcidFragment])
df_corrected_MIDs = @transform(gdfs_corrected, 
    :MID = :Response ./ sum(:Response),
    :CorrectedMID = :CorrectedResponse ./ sum(:CorrectedResponse))
select(df_corrected_MIDs, :Sample, :AminoAcidFragment, :Response, :CorrectedResponse, :MID, :CorrectedMID)
## And of course the finished DataFrame can then be saved.
## CSV.write("Corrected_MIDs.csv", df_corrected_MIDs)


# To visualize what `isotope_correction()` does we can plot the MIDs before and after correction.
using ColorSchemes, CairoMakie

gdfs = @groupby(df_corrected_MIDs, :Sample, :AminoAcidFragment)

colors = colorschemes[:bamako][1:43:256];
group_color = [PolyElement(color = color, strokecolor = :transparent)
    for color in colors]
labels = ["m+0", "m+1", "m+2", "m+3", "m+4", "m+5"]

fig = Figure();

function myplot(group, fig, colors)
    name = group.AminoAcidFragment[1]
    sample = group.Sample[1]
    ax = Axis(fig, 
        xticks = (1:2, ["uncorrected", "corrected"]),
        title = "Sample $sample $name MID")
    barplot!(ax, repeat([1, 2], inner = length(group.MID)), 
        cat(group.MID, group.CorrectedMID, dims = 1),
        stack = repeat(group.Mplus, 2),
        colormap = colors,
        color = repeat(group.Mplus .+ 1, 2))
end
myplot(gdfs[1], fig[1,1], colors[1:5])
myplot(gdfs[2], fig[1,2], colors)
myplot(gdfs[3], fig[2,1], colors[1:5])
myplot(gdfs[4], fig[2,2], colors)

Legend(fig[1:2,3], group_color, labels, framevisible = false)

fig
## save("Methionine_MIDs.png", fig)
