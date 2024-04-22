# Thermal Structure and Composition of Jupiter's Great Red Spot from JWST/MIRI
Supplementary material for Harkett et al. (2024).

Jake Harkett, Leigh N. Fletcher, Oliver R.T. King, Michael T. Roman, Henrik Melin, Heidi B. Hammel, Ricardo Hueso, Michael H. Wong, Stefanie N. Milam, Glenn S. Orton, Patrick G.J. Irwin, Imke de Pater, Thierry Fouchet, Pablo Rodríguez-Ovalle, Patrick M. Fry, & the GTO/ERS team (2024), Thermal Structure and Composition of Jupiter's Great Red Spot from JWST/MIRI, (in prep).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10785462.svg)](https://doi.org/10.5281/zenodo.10785462)

## Contents

### context_data
Contains Hubble visual context data taken on 2022-07-28, an amateur cylindrical map (Miyazaki, I.) compiled from visible observations made by; Mike Wölle (Austria), Jacques van der Meer (France) and Ed Grafton (U.S.A) on 2022-08-15 and mid-Infrared VLT/VISIR data taken on 2022-08-08.

### figures

Quick-look plots and figures used for Harkett et al. (2024) generated for the JWST/MIRI data and resulting atmospheric retrievals.

### flat_field

Contains the derived flat-field frames calculated from the observed GRS data ([Fletcher et al., 2023](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2023JE007924); [King et al., 2023](https://iopscience.iop.org/article/10.3847/2515-5172/ad045f)).

### nemesis_files

Prior data and codes used to run the NEMESIS atmospheric retrievals on the raw MIRI data.

### scripts

All other codes used in this study including those used to generate the plots and visuals stored in ```/figures/```.

### spx_files

Spectrum (SPX) files containing spectral data for each epoch, MIRI tile and band (ch1-short to ch3-medium).

## See also

The [pipelines repository](https://github.com/JWSTGiantPlanets/pipelines) ([King et al., 2023](https://iopscience.iop.org/article/10.3847/2515-5172/ad045f)) contains the reduction pipeline code used to calibrate the raw MIRI data before performing the NEMESIS atmospheric retrievals.

The [JWST Giant Planets repository](https://github.com/JWSTGiantPlanets) contains further tools for calibrating, reducing and analysing data from JWST/MIRI and JWST/NIRSpec.
