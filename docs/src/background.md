```@meta
CurrentModule = NaturalIsotopeCorrection
```

# Background

Documentation for [NaturalIsotopeCorrection](https://github.com/vm-vh/NaturalIsotopeCorrection.jl).


Metabolic flux analysis (MFA) is a powerful method used to assess intracellular metabolic activity by estimating 
the fluxes of different metabolic pathways. This approach is often more informative than simply measuring 
metabolite concentrations as it helps to understand how cells respond to different conditions. 
The stable isotope of carbon, carbon-13 (<sup>13</sup>C), is commonly used as a tracer in MFA to investigate metabolic 
fluxes, though other isotopes like nitrogen-15 (<sup>15</sup>N) can be used. These isotopes do not significantly alter 
the chemical properties of the molecules they label, making them ideal for tracking metabolic processes. 
However, naturally occurring isotopes can affect the accuracy of the isotope tracing, 
as their presence interferes with the signals obtained from the isotopically labeled substrates.

When conducting isotope labeling experiments, it is critical to distinguish between isotopes introduced into 
the system experimentally and those naturally present at the start of the experiment. The natural abundance (NA) 
of stable isotopes can falsify the results of isotope labeling studies if not corrected appropriately. [^1] 
For example, carbon naturally exists as  approximately 98.93% <sup>12</sup>C and 1.07% <sup>13</sup>C. 
Ignoring the NA of <sup>13</sup>C could lead to errors in the interpretation of metabolic fluxes, especially 
when using gas chromatography mass spectrometry (GCMS) as the derivatization necessary for gas chromatography 
leads to more incorporated atoms.

While early correction methods were based on the assumption that mass isotopomer distributions (MIDs) 
of labeled standards were simply shifted versions of the unlabeled standards, it is much more accurate 
to account for the non-linear distribution of stable isotopes, or "skew," caused by isotopic enrichment from 
experimental tracers. 
This package implements this “skewed” correction method, using multinomial probability theory to account 
for the non-linear distribution of naturally occurring isotopes. A least-squares optimization function ensures 
that fractional abundances remain positive after correction is used and corrects for noise, 
similarly to IsoCor. [^2] In addition to correcting for the natural abundances of all isotopes, correction for 
the purity of the tracer substrate is also implemented.

[^1]	Midani, F. S.; Wynn, M. L.; Schnell, S. The Importance of Accurately Correcting for the Natural Abundance 
of Stable Isotopes. Anal. Biochem. 2017, 520, 27–43. https://doi.org/10.1016/j.ab.2016.12.011.

[^2]	Millard, P.; Letisse, F.; Sokol, S.; Portais, J.-C. IsoCor: Correcting MS Data in Isotope Labeling 
Experiments. Bioinformatics 2012, 28 (9), 1294–1296. https://doi.org/10.1093/bioinformatics/bts127.


