# Quick guide

This GitHub repository provides the "IPECAD Open-Source Model framework" for cost-effectiveness analysis of interventions in early Alzheimer's disease.

The following versions are hosted with some features highlighted (e.g., R/spreadsheet/Shiny, Markov/microsimulation, analyses for publication/conference; for detailed features see version history below): 

- v2.4.0 (2025): Alzheimer Europe 2024 (conference), CTAD 2024 (conference). 
- v2.3.0 (2024): Updated framework with ICER/AD-ACE cross-validation (R, shiny, spreadsheet; publication; https://ronhandels.shinyapps.io/ipecad/). 
- v2.2.0 (2024): SveDem model replication IPECAD workshop 2023 (publication). 
- v2.1.0 (2023): SveDem model replication ISPOR 2023 (conference). 
- v2.0.1 (2023): SveDem model replication update. 
- v2.0.0 (2023): SveDem model replication. 
- v1.0.0 (2019): Markov multi-domain (spreadsheet, publication). 

This readme contains an installation guide and version history. 

Releases contain stable versions. The branch "main" is for development and fune tuning before release. The branch "develop" is for development, experimenting and exploring new features, it is not necessarily stable. Releases are stable versions often associated to specific applications of the model. 

# Installation guide R version (for inexperienced R/github users)

- Download files: go to https://github.com/ronhandels/IPECAD, then click "Code" (green button) and then "Download ZIP" on PC. 
- Unzip folder. 
- Move unzipped folder to desired location on PC. 
- Make sure to have R (https://cran.r-project.org/mirrors.html) and ideally RStudio (https://www.rstudio.com) installed. 
- Open R file named "IPECAD open-source model.R" in R (or RStudio)
- Install `dampack` by running the code `install.packages("dampack")` in R (this code is located under the heading `# MANUAL PREPARATION #`; to activate the code remove the `#` at the beginning of the line).
- Set the working directory to the file path where the unzipped folder is stored (this code is located under the heading `# MANUAL PREPARATION #`; change the code manually or alternatively in RStudio go to the menu "Session" then "Set Working Directory" then "To Source File Location" (or to "Choose Directory..." and choose the location of the unzipped folder)).
- Run the code by sourcing it (in R go to the menu "File" then "Source R code..."; in RStudio go to the menu "Code" then "Source")
- We recommend familiarizing with the description of the model using the details related to the releases described below.

# Cite this work

Name: IPECAD open-source model <https://github.com/ronhandels/IPECAD>. 

Citation: Handels et al. 2024 https://doi.org/10.1016/j.jval.2024.07.009

# Version details

## v2.4.0 (2024; released)

This version is an application of v2.3.0 for additional analyses presented at Alzheimer Europe 2024 and CTAD 2024. For additional information on Alzheimer Europe 2024 see abstract [folder "additional_documentation"] and presentation slides [folder "additional_documentation"]. For additional information on CTAD see abstract [folder "additional_documentation"] and presentation slides [folder "additional_documentation"]. General details on the open-source framework can be found in v2.3.0 (Handels et al. \[2024\] https://doi.org/10.1101/2024.04.05.24305373). 

Features: same as v2.3.0; with simplified flow in model functions. 

This version has been developed and/or supported by: 

- Ron Handels
- Anders Wimo
- Bengt Winblad
- Linus Jonsson
- Sabine Grimm
- Anders Skoldunger

## v2.3.0 (2024; released)

This version has been developed as an open-source framework with details found here: Handels et al. \[2024\] https://doi.org/10.1101/2024.04.05.24305373. In this publication we describe the new IPECAD open-source model framework (version 2) for the health-economic evaluation of early AD treatment and apply it in 3 use cases for AD lecanemab treatment: 1) cross-validating an existing model with a similar structure (ICER), 2) cross-validating an existing model with a more complex structure (AD-ACE) and 3) assessing additional uncertainty scenarios. 

Features: R; shiny; spreadsheet; cycle time; calibrate time shift. 

This version has been developed and/or supported by: 

- Ron Handels (main developer, Maastricht University - Netherlands & affiliated to Karolinska Institutet - Sweden)
- William L. Herring (RTI Health Solutions - USA & affiliated to Karolinska Institutet - Sweden)
- Sabine Grimm (Maastricht University Medical Centre - Netherlands)
- Anders Sköldunger (Karolinska Institutet - Sweden)
- Bengt Winblad (Karolinska Institutet - Sweden)
- Anders Wimo (Karolinska Institutet - Sweden)
- Linus Jönsson (Karolinska Institutet - Sweden)

## v2.2.0-beta (2024; released)

This version is used for the IPECAD workshop model cross-comparison 2023. Details can be found here: 

-  ISPOR 2023 abstract: <https://www.ispor.org/heor-resources/presentations-database/presentation/euro2023-3788/131246>
-  ISPOR 2023 poster: <https://osf.io/7xbez>
-  IPECAD workshop details: <https://ipecad.org/workshop/>
-  Publication: Handels et al. <https://doi.org/10.1016/j.jval.2024.09.006>

This version extends v2.1.0 with new features (transitions to institutionalization, table/graph outcomes), new inputs (estimate of the transition from MCI to mild dementia, estimate for state-specific costs in community and institutional setting) and fixes (half-cycle correction before discounting). 

Features: R. 

This version has been developed by: 

- Ron Handels (main developer, Maastricht University - Netherlands)
- Anders Wimo (Karolinska Institutet - Sweden)
- Linh Nguyen (Maastricht University - Netherlands)

## v2.1.0-beta (2023; released)

This version was build for the abstract submitted to ISPOR Europe 2023. For additional information see [presentation slides](https://github.com/ronhandels/IPECAD/blob/main/additional_documentation/presentation_2023_11_14_IPECAD_open_source_model_ISPOR_v4.pdf) presented on ISPOR Europe 2023 and [session details](https://www.ispor.org/conferences-education/conferences/upcoming-conferences/ispor-europe-2023/program/program/session/euro2023-3781/17212). 

Features: R; headroom; probabilistic. 

This version has been developed and/or supported by:

-   Ron Handels (main developer, Maastricht University - Netherlands)
-   Anders Wimo
-   Anders Sköldunger
-   Ashley Tate
-   Bengt Winblad
-   Daphne Silvertand
-   Linh Nguyen
-   Linus Jönsson
-   Sabine Grimm
-   Sandar Aye
-   Will Herring

## v2.0.1-alpha (2023)

We replaced any (non-rounded input) estimates or information not in the public domain with (rounded) input estimates or information available in the public domain. Specifically, this refers to the following aspects:

-   The original 'SveDem' model used a Swedish life table to apply mortality. The mortality rate was based on the weighted mean rate for men and women (weights based on the prevalence of age and sex in the Swedish population). For the model replication this was changed by using an updated Swedish life table and a sex-specific mortality rate.
-   The original 'SveDem' model multiplied the transition probability from MCI to mild dementia and from mild dementia to moderate dementia with the treatment effect expressed as a ratio. For the model replication this was changed by first transforming the transition probability to a rate, then multiplying it with the treatment effect expressed as a ratio, then transformed back to a transition probability.
-   The original 'SveDem' model used non-rounded input estimates for transition probabilities and mortality hazard ratios. For the model replication this was changed to only using estimates as publicly available in the publication (see supplemental material of earlier 'SveDem' model publication).

This version is not publicly available given its the alpha version status.

This version has been developed by:

-   Ron Handels
-   Anders Wimo
-   Linh Nguyen
-   Daphne Silvertand

## v2.0.0-alpha (2023)

An existing model \[Wimo et al. 2020: <https://www.doi.org/10.3233/JAD-191055>\] (further referred to as 'original SveDem model') was reconstructed in terms of its disease progression part. This original SveDem model was a cohort-based Markov model and consisted of health states Mild Cognitive Impairment (MCI), mild dementia, moderate dementia, severe dementia and death. Progression between states was modeled with a cycle length of 1 year simulating from age 60 to 100. See publication for a detailed acknowledgment of the developers and funders (partly funded by Merck & Co., Inc., Kenilworth, NJ, USA).

Reconstruction of the original SveDem model was done by the original developers Anders Wimo and Ron Handels using original non-rounded input estimates available in the original model format TreeAge and Excel. The reconstructed model produced the same model outcomes in terms of mean person-years per person per state and alive (table 1 and table 4 from earlier mentioned 'SveDem' model publication) when rounded to 2 decimal points, except for person-years in severe dementia after 40 years which had an absolute deviation of 0.01 person-years. Anders Wimo and Ron Handels considered this a sufficient reflection of internal validity. No specific funding was received for this reconstruction.

This version is not publicly available given the dependency on non-publicly available information and its alpha version status.

This version has been developed by:

-   Ron Handels
-   Anders Wimo

## v1 (2019)

The IPECAD Open-Source Model v1 - Multi-Domain is an open-source model described in a publication \[Green et al. 2019: <https://doi.org/10.1016/j.jalz.2019.05.004>\] and can be freely requested at <https://www.ipecad.org/open-source-model>. This v1 is a different model independent from v2.

# Acknowledgment

We acknowledge all members of the IPECAD group \[<https://www.ipecad.org>\].

We acknowledge all (future) contributors to the model.

We acknowledge the developers of the `dampack` (<https://github.com/DARTH-git> and <http://darthworkgroup.com>).
