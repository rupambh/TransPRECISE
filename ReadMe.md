# TransPRECISE: Personalized Network Modeling of Pan-cancer Patient and Cell Line Interactome

Hi! Welcome to our repository for the TransPRECISE project!

## General Comments and Guides

The directories and files are, in general, named as self-explanatorily as possible! The .R script files are organized in a way so that all the dependencies required to run them are in the same directory as well. If you are unsure of the purpose of any of the files and/or think that some necessary scripts/data files are missing, let us know at [this email](mailto:rupamb@umich.edu).

## Directories and Contents

1. The **Data** folder contains all the base data files used to generate all our results. All further used .rda files are generated from some combination of these files.
2. The **Figures** folder contains all the codes for generating the figures in our results section, organized and numbered in the same order as in the paper. For cases where parts of the figures are generated separately, we used [Omnigraffle 7](https://store.omnigroup.com/main/omnigraffle) for merging and editing.
3. The **PRECISE** folder contains the codes to fit the Bayesian graphical regression models in order to compute, for each cancer type (patients/cell lines), the cancer-specific pathway networks, and the sample-specific networks, activity statuses, and scores. For this purpose, we use the [PRECISE](https://github.com/MinJinHa/PRECISE) pipeline developed by Ha et al., 2018.
4. The **Revisions** folder contains three things - the codes used to run a confirmatory analysis using the drug exposure data available from the [Gene-Drug Interactions for Survival in Cancer (GDISC)](https://gdisc.bme.gatech.edu/) project, the codes and results for our case study on the head and neck patient and cell lines samples, and a set of tutorials to illustrate the Bayesian graphical regression algorithm.
5. The **Tables** folder contains all the codes for generating the supplementary tables, organized and numbered in the same order as in the paper. The R outputs are edited and formatted using [Microsoft Office 365](https://www.office.com/).

## R Shiny App

The R Shiny app for interactive visualizations of our results is hosted [here](http://rupamb.shinyapps.io/transprecise).

## Citation

If you are using TransPRECISE codes or the Shiny app, please cite the following.

Bhattacharyya, R., Ha, M. J., Liu, Q., Akbani, R., Liang, H. & Baladandayuthapani, V. (2019+). Proteomics-based Network Modeling of Pan-cancer Human and Cell line Interactome. *Submitted to JCO CCI, under revision*.
