"""
Describe me
"""
module NaturalIsotopeCorrection

using DocStringExtensions
using DSP
using LBFGSB

include("types.jl")
include("read.jl")
include("matrices.jl")
include("correction.jl")

export isotope_correction

end
