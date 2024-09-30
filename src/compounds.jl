"""
$(TYPEDSIGNATURES)

Outputs a tuple containing two elements:
- A dictionary of all elements in the fragment and their respective amounts.
- The number of possible isotopologues (which is equal to the number of possibly labeled carbon atoms + 1).

# Arguments
- `formula::String`: full chemical formula of the fragment to be corrected e.g. "C3H3O3LabC3" (or "C3H3O3" in conjunction with label)
- `label::String`: labeled element and amount e.g. "C3", can be ignored if the formula contains LabC
"""
function fragment(formula::String, label::String)::Tuple{Dict, Int}
    if label == ""
        split_formula = split(formula, "LabC", keepempty = false)
        length(split_formula) === 1 ? throw(ArgumentError("Only carbon is currently implemented as a tracer, please make sure your formula is a string containing LabC and a number 
                      or use the label keyword.")) : nothing
        formula = split_formula[1]
        n_isotopologues = parse(Int,  split_formula[2]) + 1
    elseif label[1] != 'C'
        throw(ArgumentError("Only carbon is currently implemented as a tracer, please make sure your label is a string containing only C and a number."))
    else
        n_isotopologues = parse(Int,  label[2:end]) + 1
    end
    
    elements = split(formula, r"[0-9]", keepempty = false)
    amounts = split(formula, r"[aA-zZ]", keepempty = false)
    fragment_dict = Dict(elements .=> map(x -> parse(Int, x), amounts))
    fragment_dict["C"] = fragment_dict["C"] - (n_isotopologues - 1)

    fragment_dict["C"] < 0 ? throw(ArgumentError("The number of labeled tracer atoms cannot be larger than the total number of tracer atoms.")) : nothing

    return fragment_dict, n_isotopologues
end

#= TODO: implement this

abstract type AbstractLabeledCompound end

struct LabeledCompound <: AbstractLabeledCompound
    LabC::Int = 0
    H::Int = 0
    C::Int = 0
    N::Int = 0
    O::Int = 0
    S::Int = 0
    Si::Int = 0
end
=#
