<!-- output: html_document -->
<!-- word_document: default -->
<!-- output:  -->

<!-- md_document: --> <!--   variant: markdown_github -->

This GitHub repository provides an open-source model for
cost-effectiveness analysis of Alzheimer’s disease interventions.

The model is available in 3 formats under a public license that allows
obtaining the source code and making changes.

| Format      | Comment                                     | Link                                                                       |
|:----------------|:---------------------------|:---------------------------|
| R           | using base R and `dampack` package          | <https://github.com/ronhandels/ipecad/> \> file `AD open-source model.R`   |
| spreadsheet | deterministic only                          | <https://github.com/ronhandels/ipecad/> \> file `AD open-source model.ods` |
| Online      | deterministic only with selection of inputs | <https://ronhandels.shinyapps.io/ipecad/>                                  |

Generally, the formats are based on the same underlying inputs list and
strategy and scenario functions (see model code for details). See github
for details on versioning and potential differences between them.

## Installation tutorial

-   go to <https://github.com/ronhandels/ipecad/> \> click on ‘\<\>
    Code’ (green button) \> click on ‘Download ZIP’ \> unzip all files
-   install R
-   install RStudio (optional but suggested)
-   install `dampack` by running the code `install.packages("dampack")`
-   open file ‘AD open-source model.R’ in RStudio (or R)
-   in the code under chapter ‘TECHNICAL PREPARATION’ change the working
    directory `~/GitHub/IPECAD` to the directory you have stored the
    file ‘AD open-source model.R’ in to make sure the life table in the
    subfolder is correctly referred to.
-   we recommend familiarizing with the general description of the model
    provided below and the description of its input sources detailed
    below.

## Background

Limited health-care resources require choices on how to distribute them.
Cost-effectiveness analysis assesses the value an intervention to inform
such choice. A decision-analytic model can be used when available
empirical health-economic data are limited, for example to extrapolate
beyond a randomized controlled trial follow-up period. Transparency and
credibility of decision-models are key for their results to be used by
decision-makers.

Recent drug developments in Alzheimer’s disease (AD) have increased
attention for decision-analytic models to assess their value. We aim to
develop an open-source AD decision-model framework for the
health-economic evaluation of AD interventions. We believe an
open-source model supports transparency and credibility, and can be used
to externally validate results generated by other decision-analytic
models in AD. We welcome attributions to the model on github.

Here we describe the methods and code of the open-source model for
cost-effectiveness analysis of Alzheimer’s disease interventions.

## Method

A cohort state-transition Markov model has been developed. Initial
developement consisted of reconstructing an existing model in terms of
its disease progression part. Then, inputs were changed and new features
were added.

### Features and assumptions

The model has the following features and assumptions:

-   The model simulates 2 strategies: standard of care (SOC) and
    intervention (INT); further referred to as a scenario (e.g., base
    case scenario).
-   The markov states include mild cognitive impairment (MCI), mild
    dementia (MIL), moderate dementia (MOD), severe dementia (SEV) and
    death (DTH); with MCI and mild dementia split into on- and
    off-treatment. Dementia health states were defined by the
    Mini-Mental State Examination (MMSE) categorized into mild (30-21),
    moderate (20-10) and severe (0-9).
-   Transitions between health states other than death are simulated
    time-independent.
-   Non-death transition are conditional on being alive. Some
    non-sequential and backtransitions are allowed. On/off treatment
    conditions are independent of health state and unidirectional (only
    on to off).
-   Mortality is implemented by multiplying an age- and sex-specific
    life table (user-defined) with the relative risk for death by
    disease severity state.
-   A cycle length of 1 year is employed with half-cycle correction
    operationalized as taking the mean of the health-economic outcomes
    at each 2 subsequent cycles.
-   The time horizon is variable, up to the age of 99.
-   The base population (user-defined) is defined by age (single age
    between 50-99) and sex (male or female), both used to select the
    corresponding age- and sex-specific mortality rate from a general
    population life table at each cycle.
-   Treatment effect is implemented as a relative risk with which
    forward transition rates (user-selected) were multiplied with. In
    the standard of care strategy all start off treatment, in the
    intervention strategy all start in on treatment.
-   Discounting (user-defined) is applied on the half-cycle corrected
    estimates.
-   The model runs for male or female separately. It is up to the
    end-user to provide weighted mean results.

### Input estimates of the base case scenario

The following input estimates were used for the base-case scenario:

| Parameter            | Value                                                           | Source                                                         | Notes                                                                                                                                                                                                                                                                                                                                |
|:------------|:------------|:------------|:----------------------------------|
| age_start            | 70                                                              | user-defined                                                   | age at start                                                                                                                                                                                                                                                                                                                         |
| age_end              | 99                                                              | user-defined                                                   | age up to which model is run (reflecting time horizon)                                                                                                                                                                                                                                                                               |
| sex                  | “female”                                                        | user-defined                                                   | sex                                                                                                                                                                                                                                                                                                                                  |
| p.mci_mil            | 0.206                                                           | \[Wimo, 2020: <https://doi.org/10.3233/jad-191055>\]           | transition probability MCI to mild dementia                                                                                                                                                                                                                                                                                          |
| p.mci_mod            | 0                                                               | idem                                                           | transition probability MCI to moderate dementia                                                                                                                                                                                                                                                                                      |
| p.mci_sev            | 0                                                               | idem                                                           | transition probability MCI to severe dementia                                                                                                                                                                                                                                                                                        |
| p.mil_mci            | 0                                                               | idem                                                           | transition probability mild dementia to MCI                                                                                                                                                                                                                                                                                          |
| p.mil_mod            | 0.293                                                           | idem                                                           | transition probability mild dementia to moderate dementia                                                                                                                                                                                                                                                                            |
| p.mil_sev            | 0.001                                                           | idem                                                           | transition probability mild dementia to severe dementia                                                                                                                                                                                                                                                                              |
| p.mod_mil            | 0.087                                                           | idem                                                           | transition probability moderate dementia to mild dementia                                                                                                                                                                                                                                                                            |
| p.mod_sev            | 0.109                                                           | idem                                                           | transition probability moderate dementia to severe dementia                                                                                                                                                                                                                                                                          |
| p.sev_mil            | 0                                                               | idem                                                           | transition probability severe dementia to mild dementia                                                                                                                                                                                                                                                                              |
| p.sev_mod            | 0.196                                                           | idem                                                           | transition probability severe dementia to moderate dementia                                                                                                                                                                                                                                                                          |
| m.r.mortality        | m.mortality_rate_US                                             | \[<https://www.cdc.gov/nchs/data/nvsr/nvsr70/nvsr70-19.pdf>\]  | general population mortality rate by age and sex (life table)                                                                                                                                                                                                                                                                        |
| hr.mort_mci          | 1                                                               | assumption                                                     | hazard ratio mortality MCI assumed same as general population                                                                                                                                                                                                                                                                        |
| hr.mort_verymilddem  | 1.82                                                            | \[Andersen, 2010: <https://doi.org/10.1159/000265553>\]        | hazard ratio mortality very mild dementia                                                                                                                                                                                                                                                                                            |
| hr.mort_mil          | 1.318                                                           | \[Wimo, 2020: <https://doi.org/10.3233/jad-191055>\]           | hazard ratio mortality mild dementia (hazard ratio by dementia state compared to very mild dementia \[Wimo, 2020: <https://doi.org/10.3233/jad-191055>\] multiplied with HR or very mild compared to no dementia \[Andersen, 2010: <https://doi.org/10.1159/000265553>\])                                                            |
| hr.mort_mod          | 2.419                                                           | idem                                                           | hazard ratio mortality moderate dementia                                                                                                                                                                                                                                                                                             |
| hr.mort_sev          | 4.267                                                           | idem                                                           | hazard ratio mortality severe dementia                                                                                                                                                                                                                                                                                               |
| rr.tx_mci_mil        | 0.75                                                            | user-defined                                                   | relative risk treatment effect transition MCI to mild dementia                                                                                                                                                                                                                                                                       |
| rr.tx_mci_mod        | 1                                                               | user-defined                                                   | relative risk treatment effect transition MCI to moderate dementia                                                                                                                                                                                                                                                                   |
| rr.tx_mci_sev        | 1                                                               | user-defined                                                   | relative risk treatment effect transition MCI to severe dementia                                                                                                                                                                                                                                                                     |
| rr.tx_mil_mod        | 0.75                                                            | user-defined                                                   | relative risk treatment effect transition mild dementia to moderate dementia                                                                                                                                                                                                                                                         |
| tx_waning            | 0.05                                                            | user-defined                                                   | annual treatment waning effect expressed as relative reduction in treatment effect                                                                                                                                                                                                                                                   |
| p.discontinuation1   | 0.1                                                             | user-defined                                                   | transition probability from on to off treatment in year 1 (discontinuation)                                                                                                                                                                                                                                                          |
| p.discontinuation_x  | 0.1                                                             | user-defined                                                   | transition probability from on to off treatment after year 1 (discontinuation)                                                                                                                                                                                                                                                       |
| tx_duration          | 7                                                               | user-defined                                                   | maximum treatment duration                                                                                                                                                                                                                                                                                                           |
| p.starting_state_mci | 1                                                               | user-defined                                                   | proportion starting population in state MCI (1 minus this proportion starts in mild dementia); in the ‘soc’ strategy all start off-treatment, in ‘int’ strategy all start on-treatment                                                                                                                                               |
| u.mci                | 0.73                                                            | \[Green, 2019: <https://doi.org/10.1016/j.jalz.2019.05.004>\]  | utility MCI                                                                                                                                                                                                                                                                                                                          |
| u.mil                | 0.69                                                            | idem                                                           | utility mild dementia                                                                                                                                                                                                                                                                                                                |
| u.mod                | 0.53                                                            | idem                                                           | utility moderate dementia                                                                                                                                                                                                                                                                                                            |
| u.sev                | 0.38                                                            | idem                                                           | utility severe dementia                                                                                                                                                                                                                                                                                                              |
| c.mci                | (1254 + 222) \* 12 \* (1-0 ) + (1254 + 8762) \* 12 \* 0         | \[Tahami, 2023: <https://doi.org/10.1007/s40120-023-00460-1>\] | costs MCI (build up as montly costs in patient health and social care by care setting (community/residential) multiplied by 12 (annual costs) \[Tahami, 2023: <https://doi.org/10.1007/s40120-023-00460-1> table 2\] and multiplied by proportion in setting \[Tahami, 2023: <https://doi.org/10.1007/s40120-023-00460-1> table 1\]) |
| c.mil                | (1471 + 410) \* 12 \* (1-0.038) + (1471 + 8762) \* 12 \* 0.038  | idem                                                           | costs mild dementia                                                                                                                                                                                                                                                                                                                  |
| c.mod                | (1958 + 653) \* 12 \* (1-0.110) + (1958 + 8762) \* 12 \* 0.110  | idem                                                           | costs moderate dementia                                                                                                                                                                                                                                                                                                              |
| c.sev                | (2250 + 1095) \* 12 \* (1-0.259) + (2250 + 8762) \* 12 \* 0.259 | idem                                                           | costs severe dementia                                                                                                                                                                                                                                                                                                                |
| c.Tx                 | 10000                                                           | user-defined                                                   | costs treatment                                                                                                                                                                                                                                                                                                                      |
| c.Tx_diagnostics1    | 2000                                                            | user-defined                                                   | costs diagnostics cycle 1 (not half-cycle corrected)                                                                                                                                                                                                                                                                                 |
| discount_QALY        | 0.035                                                           | user-defined                                                   | discount rate QALY                                                                                                                                                                                                                                                                                                                   |
| discount_COST        | 0.035                                                           | user-defined                                                   | discount rate costs                                                                                                                                                                                                                                                                                                                  |
| wtp                  | 40000                                                           | user-defined                                                   | willingness to pay                                                                                                                                                                                                                                                                                                                   |

<!-- !!!TO-DO: @ daphne/linh: please tidy up this table and copy this tables description to the excel file and also as coding notes (#) to the code below and the R file on github.  -->

Notes on input estimates for the base case scenario:

-   p.mci_mil: transition probability between MCI and dementia was
    obtained from a rate reflecting AD high likelihood group \[Vos et
    al. 2015: <https://doi.org/10.1093/brain/awv029>\]. Transitions were
    assumed to mild stage of dementia (not moderate or severe).
-   p.mil_mod, p.mil_sev… etc: transition probabilities between dementia
    states were obtained from predicted values of an ordered probit
    model fitted to data from the Swedisch Dementia Registry (SveDem)
    longitudinal cognitive status after controlling for drop-out
    \[Handels et al. 2020: <https://doi.org/10.1002/alz.12050>\].
    Cognitive status was reflected by the Mini-Mental State Examination
    (MMSE) and categorized as mild dementia (21-30), moderate dementia
    (10-20) or severe dementia (0-9) health state.

## Basic overview of model code (R version)

The model code is build-up as follows:

The model is run using 2 functions: run function `f.run_strategy` and
run function `f.run_scenario`. The second function (`f.run_scenario`)
includes a loop over all strategies by calling them one by one. Overall,
the code follows these steps:

-   A: function to run a scenario
-   B: prepare and initialize objects to store scenario and strategy
    outcomes
-   C: run each strategy in a loop
-   D: prepare inputs to be used in each strategy
-   E: run preparations specific for the intervention strategy
-   F: store newly created or updated inputs to be used for the function
    to run a single strategy
-   G: run the strategy
    -   G1: prepare transition probability matrix
    -   G2: some checks
    -   G3: initialize objects to store strategy outcomes
    -   G4: starting state
    -   G5: markov multiplication by looping over cycles
    -   G6: multiply states with utility and cost estimates
    -   G7: apply half-cycle correction
    -   G8: apply discounting
    -   G9: store outcomes to be wrapped up by the ‘run scenario’
        function
-   H: store strategy results
-   I: add strategy results to scenario outcomes

## Model code

See file `AD open-source model.R` associated with comments to supports
its readibility and interpretation.

<!-- !!!TO-DO: @ daphne/linh: i will explain during phone call, but idea is to show simple examples of the difficult parts in the model code. Among which:  -->
<!-- 1. how the cycle-dependent TP matrix looks like (basically, explain how arrays work in a simple example and refer to details online R book).  -->
<!-- 2. explain how the 2 functions work with a simple example (so function `scenario` is used to prepare the inputs, then function `strategy` runs a strategy, then `scenario` takes the results and stores it (something along these lines). And refer to the available dampack R vignette and indicate we used that as an example.  -->
<!-- 3. explain the Markov principle and data dependencies:  -->
<!-- The model is build-up with the following dependencies: -->
<!-- * starting demographics (age and sex) and time are independent.  -->
<!-- * disease state 'dementia onset' is a function of time.  -->
<!-- * disease states within dementia are a function of the history if disease state (only 1 cycle back); indirectly they are a function of time.  -->
<!-- * death is a function of demographics and disease state.  -->
<!-- * QALYs & costs are a function of disease state and death.  -->
<!-- 4. add a paragraph on pre-fixes (and refer to the coding framework by the DARTH group):  -->
<!-- Standard naming for objects is 'x.y.object_name' with: x: -->
<!-- v = vector -->
<!-- m = matrix -->
<!-- a = array -->
<!-- df = data frame -->
<!-- l = list y: -->
<!-- p = probability or proportion -->
<!-- r = rate -->
<!-- rr = relative risk -->
<!-- hr = hazard ratio -->
<!-- n = number -->
<!-- u = utility -->
<!-- c = cost -->
<!-- Other naming: -->
<!-- f.name = function -->
<!-- temp.name = temporary object -->
<!-- 5. Debugging basics: explain really quickly how the debuggint function in rstudio works (search for a video and web-page you can refer to). -->
<!-- 6. explain PSA (and refer to dampack) -->
<!-- 7. refer to a video on basic principles of markov matrix multiplication (part in the code with the %*%). Refer as much as possible to the video's on youtube "Markov modelling (theory only); Decision analytic modelling in health economics": https://youtube.com/playlist?list=PLIT7NqYN7YT3LNYOoTEFchRXX7crLfYa5  -->
<!-- 8. explain how a list works -->
<!-- 9. explain how treatment waning works: temp.waning <- (1-tx_waning)^(0:(n.cycle-1)) -->
<!-- 10. explain that '# probability of remaining in the same state' needed to be updated in the intervention strategy -->

## Version details

### v2.0.0-alpha (2023)

An existing model \[Wimo et al. 2020:
<https://www.doi.org/10.3233/JAD-191055>\] (further referred to as
‘original SveDem model’) was reconstructed in terms of its disease
progression part.

This original SveDem model was a cohort-based Markov model and consisted
of health states Mild Cognitive Impairment (MCI), mild dementia,
moderate dementia, severe dementia and death. Progression between states
was modeled with a cycle length of 1 year simulating from age 60 to 100.

Reconstruction of the original SveDem model was done with support of the
original developers Anders Wimo and Ron Handels using original
non-rounded input estimates available in the original model format
TreeAge and Excel. The reconstructed model produced the same model
outcomes in terms of mean person-years per person per state and alive
(table 1 and table 4 from earlier mentioned ‘SveDem’ model publication)
when rounded to 2 decimal points, except for person-years in severe
dementia after 40 years which had an absolute deviation of 0.01
person-years. Anders Wimo and Ron Handels considered this a sufficient
reflection of internal validity.

This version is not publicly available given the dependency on
non-publicly available information and its alpha version status.

This version has been developed by:

-   Ron Handels
-   Anders Wimo

### v2.0.1-alpha (2023)

We replaced any (non-rounded input) estimates or information not in the
public domain with (rounded) input estimates or information available in
the public domain. Specifically, this refers to the following aspects:

-   The original ‘SveDem’ model used a Swedish life table to apply
    mortality. The mortality rate was based on the weighted mean rate
    for men and women (weights based on the prevalence of age and sex in
    the Swedish population). For the model replication this was changed
    by using an updated Swedish life table and a sex-specific mortality
    rate.
-   The original ‘SveDem’ model multiplied the transition probability
    from MCI to mild dementia and from mild dementia to moderate
    dementia with the treatment effect expressed as a ratio. For the
    model replication this was changed by first transforming the
    transition probability to a rate, then multiplying it with the
    treatment effect expressed as a ratio, then transformed back to a
    transition probability.
-   The original ‘SveDem’ model used non-rounded input estimates for
    transition probabilities and mortality hazard ratios. For the model
    replication this was changed to only using estimates as publicly
    available in the publication (see supplemental material of earlier
    ‘SveDem’ model publication).

This version is not publicly available given its the alpha version
status.

This version has been developed by:

-   Ron Handels
-   Anders Wimo
-   Linh Nguyen
-   Daphne Silvertand

### v2.1.0-beta

This version was build for the abstract submitted to ISPOR Europe 2023.
This version has been developed and/or supported by:

Ron Handels Anders Wimo Anders Sköldunger Ashley Tate Bengt Winblad
Daphne Silvertand Linh Nguyen Linus Jönsson Sabine Grimm Sandar Aye Will
Herring

## Acknowledgment

We acknowledge the developers of the `dampack`
\[<https://github.com/DARTH-git> and <http://darthworkgroup.com/>\].

We acknowledge all members of the IPECAD group \[www.ipecad.org\].

We acknowledge all (future) contributors to the model.

For support please use our github repository features:
<https://github.com/ronhandels/ipecad/>.
