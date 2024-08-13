module NaturalIsotopeCorrection

export isotope_correction
using DSP # needed for the convolution
using LBFGSB # for the optimization
include("internal_functions.jl")

"""
    isotope_correction(response, formula)

Outputs the corrected LC-/GCMS response, mass isotopomer distribution (MID) and mean enrichment of a fragment.
Additionally outputs the residuum when using L-BFGS-B for optimization.

# Arguments
- `response::Vector`: the response / measurement vector obtained from LC-/GCMS, must have 
    the length l+1 if l is equal to the amount of possibly labeled atoms in the fragment
- `formula::String`: chemical formula of the fragment to be corrected.
- `tracer_purity = 1.0`: purity of the tracer, the default represents 100%.
- `optimization::Bool = true`: whether to use the L-BFGS-B optimization algorithm.
# Examples
```julia-repl
julia> pyruvate = "C3H3O3LabC3"
julia> response_vec = [3500000, 1000000, 3500000, 800000]
julia> corr_response, MID, mean_enrichment, residuum = isotope_correction(response, pyruvate; tracer_purity = 0.99)
([3.6690841286892043e6, 880069.4781954776, 3.548554810172248e6, 765218.8519553166], 
[0.41398107163956405, 0.09929783371601102, 0.4003818041674765, 0.08633929047694855], 
0.5397698283704524, 
[-1.0858374563130465e-13, -2.566428685730154e-15, 1.1059455573558807e-13, -1.6554787924343889e-13])
julia> corr_response, MID, mean_enrichment = isotope_correction(response, pyruvate; tracer_purity = 0.99, optimization = false)
([3.6690841286882004e6, 880069.4781954774, 3.548554810173271e6, 765218.8519538215], 
[0.41398107163951964, 0.09929783371602752, 0.40038180416765856, 0.08633929047679423], 
0.5397698283704317)

```
"""
function isotope_correction(response, formula; tracer_purity = 1.0, optimization::Bool = true)
    CM = fragment_CM(formula, tracer_purity = tracer_purity)
    n = length(CM[1,:])
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

        fout, corr_response = lbfgsb(λ -> cost_func(λ, (trunc_response, trunc_CM)...), x, lb=0, m=5, factr=1e7, pgtol=1e-5, iprint=-1, maxfun=15000, maxiter=15000)
        corr_MID = corr_response ./ sum(corr_response)

        residuum = [v/sum(trunc_response) for v in (trunc_response - trunc_CM * corr_response)]
        mean_enrichment = sum([i*corr_MID[i] for i in eachindex(corr_MID)])/n

        return corr_response, corr_MID, mean_enrichment, residuum

    else
        corr_response = trunc_CM \ response
        corr_MID = corr_response ./ sum(corr_response)
        mean_enrichment = sum([i*corr_MID[i] for i in eachindex(corr_MID)])/n

        return corr_response, corr_MID, mean_enrichment
    end
end

end
