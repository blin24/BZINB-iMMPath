# BZINB-iMMPath

# Introduction


# Data
* The raw ZOE2.0 Kraken2/Bracken microbiome data is available in the Carolina Digital Repository (https://cdr.lib.unc.edu/concern/data_sets/5d86p890x) in the file `CDR_ZOEmicro_20220519.zip` > `ZOE2_main_Bracken/`. The data file used in the scripts in this project was processed by Dr. Hunyong Cho (https://github.com/Hunyong/) using the `C01_data.reading.bracken2021.R` script from https://github.com/Hunyong/ZOE_metagenomics_2022/tree/main/scripts. 
* The ZOE2.0 metabolite data is available in the MetaboLights repository (https://www.ebi.ac.uk/metabolights/MTBLS2215/files).

# Software
* To install the BZINB package, please run the following line:
install_github("Hunyong/BZINB")

# Description of scripts
* `process_data.R` Filter and save data in R lists of data frames for further analysis
* `simulation_lognormal.R` Simulate pairs of vectors from the bivariate lognormal distribution and calculate correlations (figure 5) - 10 replicates as an example.
* `simulation_bzinb.R` simulate pairs of vectors from the BZINB distribution and calculate correlations (figure 6) - 10 replicates as an example.
* `fit_example.R` code to fit the BZINB model on and calculate correlation for pairs of count vectors - 10 pairs at a time. 
* `spectral.m` matlab code for spectral clustering based on affinity matrix using non-negative correlations. The spectral clustering and asymmetrical clustering packages can be downloaded at https://sites.stat.washington.edu/mmp/software.html (http://www.stat.washington.edu/spectral/code/SpectraLib.tgz and http://www.stat.washington.edu/spectral/code/SpectraLib_A.tar).



