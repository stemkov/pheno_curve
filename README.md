# pheno_curve
This repository contains all of the raw data and scripts needed to reproduce the analysis for the paper "Bee phenological distributions predicted by inferring vital rates" by Stemkovski et al.

Authors:
Michael Stemkovski (m.stemkovski@gmail.com)
Aidan Fife (aidanfife@gmail.com)
Ryan Stuart (rstuart525@gmail.com)
William D. Pearse (will.pearse@imperial.ac.uk)

Abstract:
How bees shift the timing of their seasonal activity (phenology) to track favorable conditions influences the degree to which bee foraging and flowering plant reproduction overlap. While bee phenology is known to shift due to interannual climatic variation and experimental temperature manipulation, the underlying causes of these shifts are poorly understood. Most studies of bee phenology have been phenomenological and have only examined shifts of point-estimates such as first-appearance or peak timing. Such cross-sectional measures are convenient for analysis, but foraging activity is distributed across time, and pollination interactions are better described by overlap in phenological abundance curves. Here, we make simultaneous inferences about interannual shifts in bee phenology, emergence and senescence rates, population size, and the effect of floral abundance on observed bee abundance. We do this with a model of transition rates between life stages implemented in a hierarchical Bayesian framework and parameterized with fine-scale abundance time-series of the sweat bee Halictus rubicundus at the Rocky Mountain Biological Laboratory (Colorado, USA). We find that H. rubicundus’s emergence cueing was highly sensitive to the timing of snowmelt, but that emergence rate, senescence rate, and population size did not differ greatly across years. The present approach can be used to glean information about vital rates from other datasets on bee and flower phenology, improving our understanding of pollination interactions.

You can find an earlier version of the manuscript as a dissertation chapter here: https://www.proquest.com/docview/2825092784

We encourage interested readers to contact Michael Stemkovski at m.stemkovski@gmail.com. This model can be applied to many phenological distributions beyond those of bees, but the model code is somewhat complex if you're just opening it up for the first time or haven't worked with a Bayesian model before.

Folders:

raw_data contains all raw data from field work, unmodified. Bee data are found in bee_data.csv and bee_data_2021.csv. Flower data are in flower_data.csv and flower_data_2021.csv. Data on sampling effort are in sampling_effort.csv and sampling_effort_2021.csv. plot_locations.csv and rmbl_sites.kmz contain site borders that were used to generate the map in Figure 1, but are not actually used in the analysis.

clean_data contains data used in the analysis. All of the files in this folder are generated by the data_cleaning.R script, but you can skip that script by using the provided clean data files and going straight to the analysis script.

scripts contains all the cleaning and analysis code. You should run the scripts in this order: functions.R > data_cleaning.R > analysis.R > analysis_aux.R. functions.R contains model simulation and data wrangling functions that are called from other scripts. data_cleaning.R cleans the raw data to get it ready for analysis. analysis.R fits the Bayesian model, shows diagnostic plots, and makes the figures for the manuscript. analysis_aux.R runs some small analysis to get summary statistics for the manuscript Results sections. model.stan contains the Bayesian model code that's called from analysis.R. model_fit.RDS is provided as a convenient fitted model output to save users time. If you want to save several hours and just get the posteriors, load this RDS file in to your R environment and proceed with the code in analysis.R after the model fitting line. model_for_testing.stan is a slightly modified version of model.stan used for the simulation method comparison supplemental analysis.

R packaged needed to run the the scripts:
lubridate, mgsub, data.table, rstan, bayesplot, plotrix, RColorBrewer, sf (optional, for map plotting), flextable, mgcv

Notes:

The Bayesian models (the main one in analysis.R and the supplementary analysis in analysis_aux.R) take a long time to run. Rather than running the whole script all at once, we recommend skipping the model fitting lines and instead running the commented-out readRDS lines right after to load in pre-fitted model objects if you want to save several hours.

In model.stan, the linear interpolation function is borrowed from a help forum on mc-stan.org. The discussion thread is here: https://discourse.mc-stan.org/t/linear-interpolation-and-searchsorted-in-stan/13318/5 The user who put together that function is: https://discourse.mc-stan.org/u/dmi3kno/summary

If you don't want to recreate the map from Figure 1, you can skip the code that loads in the kmz file. If the kmz file doesn't load using your version of the sf R-package, this will not affect the analysis.
