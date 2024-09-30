"""
$(TYPEDSIGNATURES)

Outputs a tuple containing two elements:
- The mass distribution vector (MDV) of a fragment, calculated by convolving the vectors of 
    natural isotope abundances n times for each element with n atoms.
- The number of possible isotopologues (which is equal to the number of possibly labeled carbon atoms + 1).

# Arguments
- `formula`: full chemical formula of the fragment to be corrected e.g. "C3H3O3LabC3" (or "C3H3O3" in conjunction with label)
- `label`: labeled element and amount e.g. "C3", can be ignored if the formula contains LabC
- `isotopes`: dictionary containing the natural isotope abundances of elements, default is isotope_NAs
"""
function calculate_MDV(formula::String, label::String, isotopes::Dict)::Tuple{Vector{Float64}, Int}
    correction_formula, n_isotopologues = fragment(formula, label)

    result = [1.0]  # mass are normalized to 1
    for (element, n) in correction_formula
        for _ in 1:n
            result = conv(result, isotopes[element])
        end
    end
    return result, n_isotopologues
end

"""
$(TYPEDSIGNATURES)

Outputs the correction matrix constructed by iterative convolution of the correction vector for each column.

# Arguments
- `formula`: full chemical formula of the fragment to be corrected e.g. "C3H3O3LabC3" (or "C3H3O3" in conjunction with label)
- `label`: labeled element and amount e.g. "C3", can be ignored if the formula contains LabC
- `isotopes`: dictionary containing the natural isotope abundances of elements, default is isotope_NAs
- `tracer_purity::Float64`: purity of the tracer, default is 1.0 (corresponding to 100%)
"""
function CM_conv(formula::String, label::String, isotopes::Dict, tracer_purity::Float64)::Matrix{Float64}
    tracer_element = "C" # turn into arg once other tracers are implemented
    idx_tracer = 1 # turn into arg once other tracers are implemented

    tracer_purity > 1.0 ? throw(ArgumentError("The tracer purity cannot exceed 1.0 (corresponding to 100%).")) : nothing
    tracer_purity = [1-tracer_purity, tracer_purity]
    correction_vector, n_isotopologues = calculate_MDV(formula, label, isotopes)

    # create correction matrix
    correction_matrix = zeros((n_isotopologues, n_isotopologues))
    # Peaks to keep after convolution
    mask = [n * idx_tracer for n in 1:n_isotopologues]

    # For each atom with the same element as the tracer
    for i in 1:n_isotopologues
        column = correction_vector
        # Correction of tracer purity
        for _ in 1:i-1
            column = conv(column, tracer_purity)
        end

        # Correct natural abundance of tracer element
        for _ in 0:(n_isotopologues-i-1)
            column = conv(column, isotopes[tracer_element])
        end

        if length(column) < maximum(mask)
            column += [0.0] * (maximum(mask) - length(column))
        end

        column = [column[j] for j in mask]
        correction_matrix[:, i] = column
    end
    return correction_matrix
end
