<div align="center">
    <img src="docs/src/assets/logo_text_dark.svg?maxAge=0" width="40%">
</div>

# IsotopeCorrection.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://vm-vh.github.io/IsotopeCorrection.jl/dev/)
[![Build Status](https://github.com/vm-vh/IsotopeCorrection.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/vm-vh/IsotopeCorrection.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/vm-vh/IsotopeCorrection.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/vm-vh/IsotopeCorrection.jl)

A package for performing natural isotope abundance correction on LC-/GCMS data of 13C labeling experiments.

Please note that currently only 13C is implemented as a tracer, and correction can only be done on low resolution LC-/GCMS data.

### Example
```julia
using IsotopeCorrection

# define the response / measurement vector 
response_vec = [3500000, 1000000, 3500000, 800000]

# define the chemical formula of the fragment, e.g. pyruvate
# include all atoms introduced through derivatization if appliable
pyruvate_formula = "C3H3O3"

# define the element and number of possibly labeled atoms due to the tracer
# i.e. do not include atoms of the tracer element introduced by derivatization
pyruvate_label = "C3"

# run th correction function, the tracer purity can be adjusted using a keyword argument
corrected_response, corrected_MID, mean_enrichment, residuum = isotope_correction(response_vec,
                                                                                  pyruvate_formula,
                                                                                  label = pyruvate_label,
                                                                                  tracer_purity = 0.99)

# alternatively, pyruvate_formula could also be defined as "C3H3O3LabC3" for compatibility with existing libraries
# in this case the label keyword should be left blank
```
