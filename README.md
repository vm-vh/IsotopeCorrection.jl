<div align="center">
    <img src="docs/src/assets/logo.svg?maxAge=0" width="80%">
</div>

# IsotopeCorrection.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://vm-vh.github.io/IsotopeCorrection.jl/dev/)
[![Build Status](https://github.com/vm-vh/IsotopeCorrection.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/vm-vh/IsotopeCorrection.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/vm-vh/IsotopeCorrection.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/vm-vh/IsotopeCorrection.jl)

A package for preforming natural isotope correction on LC-/GCMS data.

### Example
```julia
using IsotopeCorrection

# define a response / measurement vector 
response_vec = [3500000, 1000000, 3500000, 800000]

# define the chemical formula of the fragment, e.g. pyruvate
pyruvate_formula = "C3H3O3LabC3"

# run th correction function, the tracer purity can be adjusted using a keyword argument
corrected_response, corrected_MID, mean_enrichment, residuum = isotope_correction(response_vec,
                                                                                  pyruvate_formula,
                                                                                  tracer_purity = 0.99)

# the correction function can also be run without optimization
corrected_response, corrected_MID, mean_enrichment = isotope_correction(response_vec,
                                                                        pyruvate_formula,
                                                                        optimization = false)
```
