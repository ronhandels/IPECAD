
# README

## INTRODUCTION

This is a new open-source model and differes from the IPECAD open-source model released in 2019 (published on https://doi-org.mu.idm.oclc.org/10.1016/j.jalz.2019.05.004 and available from www.ipecad.org/open-source-model). This new model has been initially developed by severel members of IPECAD and is hosted on the IPECAD website (see below). 

The new model is based on a reconstruction of the disease progression part of an existing model sometimes referred to as 'SveDem' (see https://www.doi.org/10.3233/JAD-191055). See VERSION HISTORY for details. 

This new model is open-source available on github (https://github.com/ronhandels/IPECAD) and is referred to on www.ipecad.org/open-source-model. 

The build-up of this model code is inspired by the dampack package vignette (https://cran.r-project.org/web/packages/dampack/vignettes/dsa_generation.html). It uses mainly base R. Graphical representation of the results relies on a package 'dampack'. 


## MODEL DEVELOPERS

Ron Handels
Anders Wimo
Linus Jonsson
Linh Nguyen
Daphne Silvertand
IPECAD open-source model group (www.ipecad.org)


## CODING REFERENCES

Standard naming for objects is 'x.y.object_name' with: 
x:
- v = vector
- m = matrix
- a = array
- df = data frame
- l = list
y: 
- p = probability or proportion
- r = rate
- rr = relative risk
- hr = hazard ratio
- n = number
- u = utility
- c = cost

Other naming: 
- f.name = function
- temp.name = temporary object



## VERSION HISTORY

### 2.1.0

. Version build for abstract ISPOR Europe 2023
. This version has been developed and supported by: 
Ron Handels
Anders Wimo
Anders Sköldunger
Ashley Tate
Bengt Winblad
Daphne Silvertand
Linh Nguyen
Linus Jönsson
Sabine Grimm
Sandar Aye
Will Herring


### 2.0.1-beta

. Update using public domain information only. 
. Not publicly available. 


### 2.0.0-alpha

. Reconstruction of original model, see details in 'RECONSTRUCTION OF ORIGINAL MODEL'. 
. Not publicly available. 
. Details: 
First, the original model (https://www.doi.org/10.3233/JAD-191055) was reconstructed with support of the original developers and use of original non-rounded input estimates available in TreeAge and Excel. The reconstructed model produced the same model outcomes mean person-years per person per state and alive (table 1 and table 4 from publication) when rounded to 2 decimal points, except for person-years in severe dementia after 40 years which had an absolute deviation of 0.01 person-years. We considered this sufficient internal validity. 
Second, we replaced any (non-rounded input) estimates or information not in the public domain with (rounded) input estimates or information available in the public domain. Specifically, this refers to the following aspects: 
- The original model used a Swedish life table, for which the Swedish mortality rate was based on the weighted mean for men and women (weights based on sex prevalence in the Swedish population). This was changed by using an updated life table and sex-specific mortality rate. 
- The original model multiplied the transition probability from MCI to mild dementia and from mild dementia to moderate dementia with the treatment effect expressed as a ratio. This was changed by first transforming the transition probability to a rate, then multiplying it with the treatment effect expressed as a ratio, then transformed back to a transition probability. 
- The original model used non-rounded input estimates for transition probabilities and mortality hazard ratios. This was changed to only using estimates as publicly available in the publication (supplemental material). 
. This version has been developed by: 
Ron Handels
Anders Wimo
Linh Nguyen
Daphne Silvertand
