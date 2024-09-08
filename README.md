# NaturalIsotopeCorrection.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://vm-vh.github.io/NaturalIsotopeCorrection.jl/dev/)
[![Build Status](https://github.com/vm-vh/NaturalIsotopeCorrection.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/vm-vh/NaturalIsotopeCorrection.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/vm-vh/NaturalIsotopeCorrection.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/vm-vh/NaturalIsotopeCorrection.jl)

A package for preforming natural isotope correction on LC-/GCMS data.

### Example
```julia
using NaturalIsotopeCorrection.jl

# define a response / measurement vector 
response_vec = [3500000, 1000000, 3500000, 800000]

# define the chemical formula of the fragment, e.g. pyruvate
pyruvate_formula = "C3H3O3LabC3"

# run th correction function, the tracer purity can be adjusted using a keyword argument
corrected_response, corrected_MID, mean_enrichment, residuum = isotope_correction(response,
                                                                                  pyruvate_formula,
                                                                                  tracer_purity = 0.99)

# the correction function can also be run without optimization
corrected_response, corrected_MID, mean_enrichment = isotope_correction(response,
                                                                        pyruvate_formula,
                                                                        optimization = false)
```
