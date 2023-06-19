
######################################## README ########################################

#################### INTRODUCTION

This is a new modeling framework as compared to IPECAD open-source model version 1.x. 
The model is based on a reconstruction of the disease progression part of an existing model sometimes referred to as 'SveDem' and described here: https://www.doi.org/10.3233/JAD-191055. 

The reconstructed model produced the same model outcomes mean person-years per person per state and alive (table 1 and table 4 from publication), except for a relatively small deviation due to rounding (largest absolute deviation was 0.006; relative deviation 0.6%). We considered this sufficient internal validity. The original model was programmed in TreeAge and exported to Excel. Model inputs were based on various pre-model data-analyses performed in Excel. The following aspects of the reconstructed model were changed to ensure the use of publicly available input estimates: 
- The original model used a Swedish life table, for which the Swedish mortality rate was based on the weighted mean for men and women (weights based on sex prevalence in the Swedish population). This was changed by using an updated life table and sex-specific mortality rate. 
- The original model multiplied the transition probability from MCI to mild dementia and from mild dementia to moderate dementia with the treatment effect expressed as a ratio. This was changed by first transforming the transition probability to a rate, then multiplying it with the treatment effect expressed as a ratio, then transformed back to a transition probability. 
- The original model used non-rounded input estimates for transition probabilities and mortality hazard ratios. This was changed to only using estimates as publicly available in the publication. 

The version is open-source available on github: https://github.com/ronhandels/IPECAD. 
The build-up of the model code is inspired by the dampack package vignette: https://cran.r-project.org/web/packages/dampack/vignettes/dsa_generation.html. It uses mainly base R. Graphical representation of the results relies on a package 'dampack'. 

For the next version we plan the following features: 
- Detailed description (R markdown). 
- Excel version.
- Loop over subgroups (to reflect heterogeneity).



#################### CODING REFERENCES

standard naming for objects is 'x.y.object_name' with
x:
- v = vector
- m = matrix
- a = array
- df = data frame
- l = list
- b = boolean
y: 
- p = probability or proportion
- u = utility
- c = cost
- r = rate
- rr = relative risk
- hr = hazard ratio
- n = number

other naming: 
- f.name # function
- temp.var # temporary variable
- var.lag # variable containing previous value
- var.0 # variable containing baseline value



#################### VERSIONING

Versioning is mainly done as major.minor.patch. 



#################### VERSION HISTORY

2.0.1

# Miscellaneous

. See introduction


# Bug fixes

. None



2.0

# Miscellaneous

. This a beta version for a new modeling framework, which is a different model from version 1.x.
. See introduction for details. 


# Bug fixes

. not applicable