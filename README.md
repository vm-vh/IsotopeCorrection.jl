# NaturalIsotopeCorrection.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://vm-vh.github.io/NaturalIsotopeCorrection.jl/stable/)
[![Build Status](https://github.com/vm-vh/NaturalIsotopeCorrection.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/vm-vh/NaturalIsotopeCorrection.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/vm-vh/NaturalIsotopeCorrection.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/vm-vh/NaturalIsotopeCorrection.jl)

A package for preforming natural isotope correction on LC-/GCMS data.

### Example
```julia
using NaturalIsotopeCorrection.jl

# define a response / measurement vector 
response_vec = [3595311, 1085606, 3520466, 808899]

# define the chemical formula of the fragment, e.g. pyruvate
chem_formula = "C3H3O3LabC3"

# run the correction function
corrMID = corrected_MID(response_vec, chem_formula, tracer_purity = 0.99)
```
The output will be a mass isotopomer distribution, corrected for the natural isotope abundances of the fragment's atoms.

```julia-repl
julia> corrMID
4-element Vector{Float64}:
 0.41527768234101536
 0.10640438947548786
 0.3930767730748226
 0.08524115510867415
```
