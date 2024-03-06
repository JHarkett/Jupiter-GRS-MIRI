# Thermal Structure and Composition of Jupiter's Great Red Spot from JWST/MIRI
Supplementary material for Harkett et al. (2024).

Jake Harkett, Leigh N. Fletcher, Oliver R.T. King, Michael T. Roman, Henrik Melin, Heidi B. Hammel, Ricardo Hueso, Michael H. Wong, Stefanie N. Milam, Glenn S. Orton, Patrick G.J. Irwin, Imke de Pater, Thierry Fouchet, Pablo Rodríguez-Ovalle, Patrick M. Fry, & the GTO/ERS team (2024), Thermal Structure and Composition of Jupiter's Great Red Spot from JWST/MIRI, (in prep).

## Contents

### context_data
Contains Hubble visual context data taken on 2022-07-28, an amateur cylindrical map (Miyazaki, I.) compiled from vivible observations made by; Mike Wölle (Austria), Jacques van der Meer (France) and Ed Grafton (U.S.A) on 2022-08-15 and mid-Infrared VLT/VISIR data taken on 2022-08-08.

### data

Contains both the raw JWST/MIRI data and the atmospheric retrieval results, there are two epochs of data; July and August. For each epoch there are three sets of observations (tiles); East, Centre and West (Note July does not contain west data). A 3-stage NEMESIS atmospheric retrieval process was used containing; stage 1, stage 2 and stage 2 Auxiliary (AUX).

```/[epoch]/[tile]/``` is the location of the raw data

```/[epoch]/[tile]/retrieval_data/``` is the location of the retrieval results for each tile

```/[epoch]/mosaics_stage[number]/``` is the location of the combined retrieval maps of temperature, ammonia, phosphine and aerosol distributions. [number] is either 1, 2 or 2_aux.

### figures

Quick-look plots and figures used for Harkett et al. (2024) generated for the JWST/MIRI data and resulting atmospheric retrievals.

### flat_field

Contains the derived flat-field frames calculated from the observed GRS data ([Fletcher et al., 2023](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2023JE007924); [King et al., 2023](https://iopscience.iop.org/article/10.3847/2515-5172/ad045f)).

### nemesis_files

Prior data and codes used to run the NEMESIS atmospheric retrievals on the raw MIRI data.

### scripts

All other codes used in this study including those used to generate the plots and visuals stored in ```/figures/```.

## See also

The [pipelines repository](https://github.com/JWSTGiantPlanets/pipelines) ([King et al., 2023](https://iopscience.iop.org/article/10.3847/2515-5172/ad045f)) contains the reduction pipeline code used to calibrate the raw MIRI data before performing the NEMESIS atmospheric retrievals.
