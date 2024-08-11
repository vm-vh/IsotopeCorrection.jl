module NaturalIsotopeCorrection

export corrected_MID
using DSP # needed for the convolution
include("internal_functions.jl")

"""
    corrected_MID(response, formula)

Outputs the corrected mass isotopomer distribution (MID) of a fragment.

# Arguments
- `response::Vector`: the response / measurement vector obtained from LC-/GCMS, must have 
    the length l+1 if l is equal to the amount of possibly labeled atoms in the fragment
- `formula::String`: chemical formula of the fragment to be corrected.
- `tracer_purity`: purity of the tracer, default is 1 (100%).
# Examples
```julia-repl
julia> response_vec = [3595311, 1085606, 3520466, 808899]
julia> corrected_MID(response_vec, "C3H3O3LabC3", tracer_purity = 0.99)
4-element Vector{Float64}:
 0.41527768234101536
 0.10640438947548786
 0.3930767730748226
 0.08524115510867415
```
"""
function corrected_MID(response::Vector, formula::String, tracer_purity = 1.0)
    full_fragment_CM = fragment_CM(formula, tracer_purity)
    corrected_response = full_fragment_CM \ response
    return corrected_response ./ sum(corrected_response)
end

end
