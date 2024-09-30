"""
$(TYPEDSIGNATURES)

Outputs the corrected response measurements, mass isotopomer distribution (MID) and mean enrichment of a fragment.
Additionally outputs the residuum when using L-BFGS-B for optimization.
Only to be used with low resolution GC-/LCMS data.

# Arguments
- `response::Vector`: the response / measurement vector obtained from LC-/GCMS, must have 
    the length n+1 if n is equal to the amount of possibly labeled atoms in the fragment
- `formula::String`: chemical formula of the fragment to be corrected.
- `label::String = ""`: labeled element and amount e.g. "C3", can be ignored if the formula contains LabC
- `isotopes::Dict = isotope_NAs`: dictionary containing the natural isotope abundances of elements
- `tracer_purity::Number = 1.0`: purity of the tracer, the default represents 100%.
- `optimization::Bool = true`: whether to use the L-BFGS-B optimization algorithm.

# Examples
```julia-repl
julia> response_vec = [3500000, 1000000, 3500000, 800000]
julia> pyruvate = "C3H3O3"
julia> pyruvate_label = "C3"
julia> corr_response, corr_MID, mean_enrichment, residuum = isotope_correction(response_vec, 
                                                                               pyruvate; 
                                                                               label = pyruvate_label)
([3.6336363843460013e6, 844166.0884382844, 3.5694537799896887e6, 780365.311750165], 
[0.4116212229745574, 0.09562780668246623, 0.4043505664463884, 0.08840040389658813], 
0.5423825378162519, 
[-1.0054050521417098e-15, -2.645802768793973e-16, -8.466568860140713e-16, -1.852061938155781e-16])

julia> pyruvate = "C3H3O3LabC3"
julia> corr_response, corr_MID, mean_enrichment, residuum = isotope_correction(response_vec, 
                                                                               pyruvate; 
                                                                               tracer_purity = 0.99)
([3.6336363843460013e6, 844166.0884382844, 3.5694537799896887e6, 780365.311750165], 
[0.4116212229745574, 0.09562780668246623, 0.4043505664463884, 0.08840040389658813], 
0.5423825378162519, 
[-1.0054050521417098e-15, -2.645802768793973e-16, -8.466568860140713e-16, -1.852061938155781e-16])
```
"""
function isotope_correction(response::Vector, formula::String; label::String = "", isotopes::Dict = isotope_NAs, tracer_purity::Number = 1.0, optimization::Bool = true)::Tuple
    CM = CM_conv(formula, label, isotopes, tracer_purity)
    n = length(CM[1,:])
    # length(response) != n ? throw(ArgumentError("The length of the response vector must equal the number of possible isotopologues 
    #                   i.e. the number of possibly labeled carbon atoms + 1.")) : nothing
    trunc_response = response[1:(n > end ? end : n)]
    trunc_CM = CM[1:length(trunc_response), 1:length(trunc_response)]

    if optimization == true
        # cost function used for optimization
        function cost_func(mid, response, CM)
            x = response - (CM * mid)
            # outputs the sum of square differences and gradient
            return (sum(x .* x), (CM' * x)*-2)
        end
        
        # the dimension of the problem
        x = fill(Cdouble(0e0), length(trunc_response))

        _, corr_response = lbfgsb(λ -> cost_func(λ, (trunc_response, trunc_CM)...), x, lb=0, m=5, factr=1e7, pgtol=1e-5, iprint=-1, maxfun=15000, maxiter=15000)
        corr_MID = corr_response ./ sum(corr_response)

        residuum = [v/sum(trunc_response) for v in (trunc_response - trunc_CM * corr_response)]
        sum(abs.(residuum)) > 1.0 ? (@warn "The residuum of the optimization is unusually high.") : nothing

        mean_enrichment = sum([i*corr_MID[i] for i in eachindex(corr_MID)])/n

        return corr_response, corr_MID, mean_enrichment, residuum

    else
        corr_response = trunc_CM \ response
        corr_MID = corr_response ./ sum(corr_response)
        mean_enrichment = sum([i*corr_MID[i] for i in eachindex(corr_MID)])/n

        return corr_response, corr_MID, mean_enrichment
    end
end
