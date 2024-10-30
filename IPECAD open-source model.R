
######################################## INFORMATION ########################################

# see readme.md on https://github.com/ronhandels/IPECAD for details



######################################## MANUAL PREPARATION ########################################

# install.packages("dampack") # remove '#' at beginning of the line and run once to install this package
setwd("~/GitHub/IPECAD") # if needed, change to the directory to the folder in which the R code and the life table folder is located



######################################## TECHNICAL PREPARATION ########################################

cat("\014") # clear console
rm(list = ls()) # clear environment
library(dampack) # load package



######################################## 1. INPUTS ########################################

######################################## 1.1. ESTIMATED ########################################

# U.S. general population life table 2019 from ssa.gov
## import life table and select only men/women and drop age column (make sure age corresponds to row number, i.e., start with age = 1)
m.lifetable_US_2019 <- as.matrix(read.csv(file="life_tables/lifetable_US_2019_ssa.csv", header=TRUE))[,c("male","female")]
## convert probability to rate
m.mortality_rate_US_2019 <- -log(1-(m.lifetable_US_2019))
## weight rate for male and female
m.mortality_rate_US_2019 <- cbind(m.mortality_rate_US_2019, weighted=NA)
m.mortality_rate_US_2019[,"weighted"] <- m.mortality_rate_US_2019[,"male"] * 0.48 + m.mortality_rate_US_2019[,"female"] * 0.52

# U.S. general population life table 2016 from cdc.gov
m.lifetable_US_2016 <- as.matrix(read.csv(file="life_tables/lifetable_US_2016.csv", header=TRUE))[,c("male","female","total")]
## convert probability to rate
m.mortality_rate_US_2016 <- -log(1-(m.lifetable_US_2016))
## weight rate for male and female
m.mortality_rate_US_2016 <- cbind(m.mortality_rate_US_2016, weighted=NA)
m.mortality_rate_US_2016[,"weighted"] <- m.mortality_rate_US_2016[,"male"] * (1-0.446) + m.mortality_rate_US_2016[,"female"] * 0.446



######################################## 1.2. MODEL INPUTS LIST ########################################

######################################## 1.2.1. INPUTS: CROSS-VALIDATION ICER ########################################

# input parameters (see readme for details)
l.inputs_icer <- list(
  v.names_state = c("mcion_c","mciof_c","milon_c","milof_c","mod_c","sev_c","mci_i","mil_i","mod_i","sev_i","dth"), # disease states: mci = mild cognitive impairment; mil = mild dementia; mod = moderate dementia; sev = severe dementia; dth = dead; x_i = living in institutional setting (without '_i' = living in community)
  v.names_strat = c("soc","int"), # strategies: soc = standard of care; int = intervention
  age_start = 71, 
  sex = "weighted", 
  p.starting_state_mci = 0.55, 
  n.cycle = 29, 
  p.mci_mil = 0.23, 
  p.mci_mod = 0, 
  p.mci_sev = 0, 
  p.mil_mci = 0.03, 
  p.mil_mod = 0.35, 
  p.mil_sev = 0.04, 
  p.mod_mil = 0.03, 
  p.mod_sev = 0.42, 
  p.sev_mil = 0, 
  p.sev_mod = 0.02, 
  p.mci_i = 0.024, 
  p.mil_i = 0.038, 
  p.mod_i = 0.110, 
  p.sev_i = 0.259, 
  m.r.mortality = m.mortality_rate_US_2019, 
  hr.mort_mci = 1.82, 
  hr.mort_mil = 2.92, 
  hr.mort_mod = 3.85, 
  hr.mort_sev = 9.52, 
  rr.tx_mci_mil = 0.69, 
  rr.tx_mci_mod = 1, 
  rr.tx_mci_sev = 1, 
  rr.tx_mil_mod = 0.69, 
  rr.tx_mil_sev = 0.69, 
  rr.tx_mci_mil_dis = 1, 
  rr.tx_mci_mod_dis = 1, 
  rr.tx_mci_sev_dis = 1, 
  rr.tx_mil_mod_dis = 1, 
  rr.tx_mil_sev_dis = 1, 
  p.tx_discontinuation1 = 0.069, 
  p.tx_discontinuation2 = 0, 
  tx_discontinuation2_begin = 2, 
  tx_duration = 29, 
  tx_waning = 0, 
  tx_waning_dis = 0, 
  u.mci_pt = 0.851 - 0.17, 
  u.mil_pt = 0.851 - 0.22, 
  u.mod_pt = 0.851 - 0.36, 
  u.sev_pt = 0.851 - 0.53, 
  u.mci_pt_i = 0.851 - 0.17, 
  u.mil_pt_i = 0.851 - 0.19, 
  u.mod_pt_i = 0.851 - 0.42, 
  u.sev_pt_i = 0.851 - 0.59, 
  u.mci_ic = -0.03, 
  u.mil_ic = -0.05, 
  u.mod_ic = -0.08, 
  u.sev_ic = -0.10, 
  u.mci_ic_i = -0.03, 
  u.mil_ic_i = -0.05, 
  u.mod_ic_i = -0.08, 
  u.sev_ic_i = -0.10, 
  u.Tx_start = -0.14 * (12/52) * 0.035, # disutility symptomatic ARIA multiplied by average duration (12 weeks) multiplied by prevalence symptomatic ARIA (3.5%)
  c.mci_hc = 6042*1.12 +  460, 
  c.mil_hc = 6042*1.56 +  965 + 0.21*365*0.333, # patient medical + informal carer medical + ChEI
  c.mod_hc = 6042*1.93 + 1544 + 0.66*365*0.333, 
  c.sev_hc = 6042*1.93 + 1930, 
  c.mci_hc_i = 6042*1.12 +  460, 
  c.mil_hc_i = 6042*1.56 +  965 + 0.21*365*0.333, 
  c.mod_hc_i = 6042*1.93 + 1544 + 0.66*365*0.333, 
  c.sev_hc_i = 6042*1.93 + 1930, 
  c.mci_sc = 0, 
  c.mil_sc = 0, 
  c.mod_sc = 0, 
  c.sev_sc = 0, 
  c.mci_sc_i = 7394*12, 
  c.mil_sc_i = 7394*12, 
  c.mod_sc_i = 7394*12, 
  c.sev_sc_i = 7394*12, 
  c.mci_ic =  69*12*32.46 + 0.204*0.049*20*52*32.46, # informal care + patient productivity loss
  c.mil_ic = 113*12*32.46 + 0.112*0.086*20*52*32.46, 
  c.mod_ic = 169*12*32.46, 
  c.sev_ic = 298*12*32.46, 
  c.mci_ic_i =  69*12*32.46*0.44 + 0.204*0.049*20*52*32.46, 
  c.mil_ic_i = 113*12*32.46*0.44 + 0.112*0.086*20*52*32.46, 
  c.mod_ic_i = 169*12*32.46*0.44, 
  c.sev_ic_i = 298*12*32.46*0.44, 
  c.Tx = 26500 + (52/2)*78.35, # drug annual wholesale acquisition cost + treatment administration frequency * administration cost
  c.Tx_start = 261.10*4 + 261.10*3*0.215, # mri cost * 3-month monitoring in year 1 + mri cost * 3 times * proportion aria
  discount_EFFECT = 0.03, 
  discount_QALY = 0.03, 
  discount_COST = 0.03, 
  wtp = 100000, 
  half_cycle_correction = TRUE
)



######################################## 1.2.1. INPUTS: CROSS-VALIDATION AD-ACE ########################################

# input parameters (see readme for details)
l.inputs_adace <- list(
  v.names_state = c("mcion_c","mciof_c","milon_c","milof_c","mod_c","sev_c","mci_i","mil_i","mod_i","sev_i","dth"), # disease states: mci = mild cognitive impairment; mil = mild dementia; mod = moderate dementia; sev = severe dementia; dth = dead; x_i = living in institutional setting (without '_i' = living in community)
  v.names_strat = c("soc","int"), 
  age_start = 73, 
  sex = "weighted", 
  p.starting_state_mci = 0.781, 
  n.cycle = 27, 
  p.mci_mil = 0.23, 
  p.mci_mod = 0, 
  p.mci_sev = 0, 
  p.mil_mci = 0.03, 
  p.mil_mod = 0.35, 
  p.mil_sev = 0.04, 
  p.mod_mil = 0.03, 
  p.mod_sev = 0.42, 
  p.sev_mil = 0, 
  p.sev_mod = 0.02, 
  p.mci_i = 0, 
  p.mil_i = 0.038, 
  p.mod_i = 0.110, 
  p.sev_i = 0.259, 
  m.r.mortality = m.mortality_rate_US_2016, 
  hr.mort_mci = 1, 
  hr.mort_mil = 2.92, 
  hr.mort_mod = 3.85, 
  hr.mort_sev = 9.52, 
  rr.tx_mci_mil = 0.69, 
  rr.tx_mci_mod = 1, 
  rr.tx_mci_sev = 1, 
  rr.tx_mil_mod = 0.69, 
  rr.tx_mil_sev = 0.69, 
  rr.tx_mci_mil_dis = 1, 
  rr.tx_mci_mod_dis = 1, 
  rr.tx_mci_sev_dis = 1, 
  rr.tx_mil_mod_dis = 1, 
  rr.tx_mil_sev_dis = 1, 
  p.tx_discontinuation1 = 0.13, 
  p.tx_discontinuation2 = 0, 
  tx_discontinuation2_begin = 27, 
  tx_duration = 27, 
  tx_waning = 0, 
  tx_waning_dis = 0, 
  u.mci_pt = 0.80, 
  u.mil_pt = 0.74, 
  u.mod_pt = 0.59, 
  u.sev_pt = 0.36, 
  u.mci_pt_i = 0.80, 
  u.mil_pt_i = 0.74, 
  u.mod_pt_i = 0.59, 
  u.sev_pt_i = 0.36, 
  u.mci_ic = -0, 
  u.mil_ic = -0.036, 
  u.mod_ic = -0.070, 
  u.sev_ic = -0.086, 
  u.mci_ic_i = -0, 
  u.mil_ic_i = -0.036, 
  u.mod_ic_i = -0.070, 
  u.sev_ic_i = -0.086, 
  u.Tx_start = -0.14 * (12/52) * 0.126, # disutility ARIA-E multiplied by average duration (12 weeks) multiplied by prevalence ARIA-E (22%)
  c.mci_hc = 1254*12, 
  c.mil_hc = 1471*12, 
  c.mod_hc = 1958*12, 
  c.sev_hc = 2250*12, 
  c.mci_hc_i = 1254*12, 
  c.mil_hc_i = 1471*12, 
  c.mod_hc_i = 1958*12, 
  c.sev_hc_i = 2250*12, 
  c.mci_sc =  222*12, 
  c.mil_sc =  410*12, 
  c.mod_sc =  653*12, 
  c.sev_sc = 1095*12, 
  c.mci_sc_i = 8762*12, 
  c.mil_sc_i = 8762*12, 
  c.mod_sc_i = 8762*12, 
  c.sev_sc_i = 8762*12, 
  c.mci_ic = 754*12 +  988*12, # caregiver healthcare + caregiver informal care
  c.mil_ic = 781*12 + 2184*12, 
  c.mod_ic = 799*12 + 3227*12, 
  c.sev_ic = 811*12 + 5402*12, 
  c.mci_ic_i = 754*12 +  435*12, 
  c.mil_ic_i = 781*12 +  961*12, 
  c.mod_ic_i = 799*12 + 1420*12, 
  c.sev_ic_i = 811*12 + 2377*12, 
  c.Tx = 0, 
  c.Tx_start = 212.14 * 5 + 0.126 * 0.78 * 212.14 * 2 + 0.126 * 0.22 * 0.91 * 796.80 + 0.126 * 0.22 * (1-0.91) * 1098.27, # monitoring (5x MRI in year 1) and ARIA-E (12.6%) being asymptomatic (78%) (2x MRI), symptomatic (22%) mild/moderate (91%) or symptomatic severe (9%)
  discount_EFFECT = 0.03, 
  discount_QALY = 0.03, 
  discount_COST = 0.03, 
  wtp = 100000, 
  half_cycle_correction = TRUE
)



######################################## 2. RUN MODEL ########################################

# The model is run using 2 functions: run function `f.run_strategy` and run function `f.run_scenario`. 
# The second function (`f.run_scenario`) includes a loop over all strategies by calling them one by one. 

# Overall, the code follows these steps: 
# A: function to run a scenario
# A1: prepare and initialize objects to store scenario and strategy outcomes
# A2: run each strategy in a loop
#   B: function to run a strategy
#   B1: inputs validity checks
#   B2: prepare inputs to be used in each strategy
#   B3: run preparations specific for the intervention strategy
#   B4: prepare transition probability matrix
#   B5: some checks
#   B6: initialize objects to store strategy outcomes
#   B7: starting state
#   B8: markov multiplication by looping over cycles
#   B9: multiply states with utility and cost estimates
#   B10: half-cycle correction to all outcomes
#   B11: discount all outcomes
#   B12: store outcomes to be wrapped up by the 'run scenario' function
# A3: store strategy results
# A4: add strategy results to scenario outcomes

# Standard naming for objects is ‘x.object_name’, with x:
# v = vector
# m = matrix
# a = array
# l = list
# df = data frame
# f = function
# temp = temporary object
# p = probability or proportion
# r = rate
# rr = relative risk or relative rate
# hr = hazard ratio
# n = number
# u = utility
# c = cost
# This is inspired by recommendations by https://github.com/DARTH-git/darthpack (https://doi.org/10.1007/s40273-019-00837-x table 3). 



######################################## 2.1. RUN STRATEGY (STEP B1-B12) ########################################

# run strategy (STEP B: function to run a strategy)
f.run_strategy <- function(l.inputs, strat) {
  with(as.list(l.inputs), { # the 'with' functions enables to call the items from the list without having to refer to the list each time (one can use 'age_start' instead l.inputs[["age_start"]])
    
    # some validity checks on input estimates (STEP B1: inputs validity checks)
    # [to be developed]
    
    # number of states
    n.state <- length(v.names_state)
    
    # convert time-independent transitions to vector of transitions (STEP B2: prepare inputs to be used in each strategy)
    v.p.mci_mil <- rep(p.mci_mil, n.cycle)
    v.p.mci_mod <- rep(p.mci_mod, n.cycle)
    v.p.mci_sev <- rep(p.mci_sev, n.cycle)
    v.p.mil_mci <- rep(p.mil_mci, n.cycle)
    v.p.mil_mod <- rep(p.mil_mod, n.cycle)
    v.p.mil_sev <- rep(p.mil_sev, n.cycle)
    v.p.mod_mil <- rep(p.mod_mil, n.cycle)
    v.p.mod_sev <- rep(p.mod_sev, n.cycle)
    v.p.sev_mil <- rep(p.sev_mil, n.cycle)
    v.p.sev_mod <- rep(p.sev_mod, n.cycle)
    
    # probability of remaining in the same state (calculated as 1 minus transitions to other states; conditional on survival)
    v.p.mci_mci <- 1 - v.p.mci_mil - v.p.mci_mod - v.p.mci_sev
    v.p.mil_mil <- 1 - v.p.mil_mci - v.p.mil_mod - v.p.mil_sev
    v.p.mod_mod <- 1 - v.p.mod_mil - v.p.mod_sev
    v.p.sev_sev <- 1 - v.p.sev_mil - v.p.sev_mod
    
    # copy for on treatment
    v.p.mcion_mci <- v.p.mci_mci
    v.p.mcion_mil <- v.p.mci_mil
    v.p.mcion_mod <- v.p.mci_mod
    v.p.mcion_sev <- v.p.mci_sev
    v.p.milon_mci <- v.p.mil_mci
    v.p.milon_mil <- v.p.mil_mil
    v.p.milon_mod <- v.p.mil_mod
    v.p.milon_sev <- v.p.mil_sev
    
    # discontinuation
    v.p.discontinuation <- rep(x=0, times=n.cycle) # initialize vector
    v.p.discontinuation[1:(tx_discontinuation2_begin-1)] <- p.tx_discontinuation1 # discontinuation up to discontinuation2 start cycle
    v.p.discontinuation[tx_discontinuation2_begin:n.cycle] <- p.tx_discontinuation2 # discontinuation at discontinuation2 start cycle onward
    v.p.discontinuation[tx_duration:n.cycle] <- 1 # maximum treatment duration implemented as discontinuation
    
    # death (subset mortality table to obtain age- and sex-specific mortality)
    v.r.dth <- m.r.mortality[age_start:(age_start+n.cycle-1), sex]
    
    # starting states
    m.trace1 <- matrix(data=0, nrow=1, ncol=n.state, dimnames=list(NULL,v.names_state))
    m.trace1[,"mciof_c"] <- p.starting_state_mci
    m.trace1[,"milof_c"] <- 1-p.starting_state_mci
    
    # strategy-specific inputs (STEP B3: run preparations specific for the intervention strategy)
    if(strat=="int") {
      
      # waning
      temp.waning <- (1-tx_waning)^(0:(n.cycle-1))
      temp.rr.tx_mci_mil <- rr.tx_mci_mil^temp.waning
      temp.rr.tx_mci_mod <- rr.tx_mci_mod^temp.waning
      temp.rr.tx_mci_sev <- rr.tx_mci_sev^temp.waning
      temp.rr.tx_mil_mod <- rr.tx_mil_mod^temp.waning
      temp.rr.tx_mil_sev <- rr.tx_mil_sev^temp.waning
      temp.waning_dis <- (1-tx_waning_dis)^(0:(n.cycle-1))
      temp.rr.tx_mci_mil_dis <- rr.tx_mci_mil_dis^temp.waning_dis
      temp.rr.tx_mci_mod_dis <- rr.tx_mci_mod_dis^temp.waning_dis
      temp.rr.tx_mci_sev_dis <- rr.tx_mci_sev_dis^temp.waning_dis
      temp.rr.tx_mil_mod_dis <- rr.tx_mil_mod_dis^temp.waning_dis
      temp.rr.tx_mil_sev_dis <- rr.tx_mil_sev_dis^temp.waning_dis
      
      # update transition probabilities treatment effect: during treatment
      v.p.mcion_mil <- 1-exp(-(-log(1-p.mci_mil) * temp.rr.tx_mci_mil)) # convert probability to rate, then multiply with treatment relative risk, then convert to probability
      v.p.mcion_mod <- 1-exp(-(-log(1-p.mci_mod) * temp.rr.tx_mci_mod)) # idem
      v.p.mcion_sev <- 1-exp(-(-log(1-p.mci_sev) * temp.rr.tx_mci_sev)) # idem
      v.p.milon_mod <- 1-exp(-(-log(1-p.mil_mod) * temp.rr.tx_mil_mod)) # idem
      v.p.milon_sev <- 1-exp(-(-log(1-p.mil_sev) * temp.rr.tx_mil_sev)) # idem
      # update transition probabilities treatment effect: after discontinuation
      v.p.mci_mil <- 1-exp(-(-log(1-p.mci_mil) * temp.rr.tx_mci_mil_dis)) # convert probability to rate, then multiply with treatment relative risk, then convert to probability
      v.p.mci_mod <- 1-exp(-(-log(1-p.mci_mod) * temp.rr.tx_mci_mod_dis)) # idem
      v.p.mci_sev <- 1-exp(-(-log(1-p.mci_sev) * temp.rr.tx_mci_sev_dis)) # idem
      v.p.mil_mod <- 1-exp(-(-log(1-p.mil_mod) * temp.rr.tx_mil_mod_dis)) # idem
      v.p.mil_sev <- 1-exp(-(-log(1-p.mil_sev) * temp.rr.tx_mil_sev_dis)) # idem
      
      # update transition probabilities of remaining in the same state
      v.p.mcion_mci <- 1 - v.p.mcion_mil - v.p.mcion_mod - v.p.mcion_sev
      v.p.milon_mil <- 1 - v.p.milon_mci - v.p.milon_mod - v.p.milon_sev
      v.p.mci_mci <- 1 - v.p.mci_mil - v.p.mci_mod - v.p.mci_sev
      v.p.mil_mil <- 1 - v.p.mil_mci - v.p.mil_mod - v.p.mil_sev
      
      # starting states
      m.trace1 <- matrix(data=0, nrow=1, ncol=n.state, dimnames=list(NULL,v.names_state))
      m.trace1[,"mcion_c"] <- p.starting_state_mci
      m.trace1[,"milon_c"] <- 1-p.starting_state_mci
      
    }
    
    # initialize time-dependent TP matrix (STEP B4: prepare transition probability matrix)
    a.TP <- array(data = 0, dim = c(n.state, n.state, n.cycle), dimnames = list(v.names_state,v.names_state,NULL))
    
    # TP matrix state: to death
    a.TP["mcion_c","dth",] <- 1-exp(-(v.r.dth * hr.mort_mci))
    a.TP["mciof_c","dth",] <- 1-exp(-(v.r.dth * hr.mort_mci))
    a.TP["milon_c","dth",] <- 1-exp(-(v.r.dth * hr.mort_mil))
    a.TP["milof_c","dth",] <- 1-exp(-(v.r.dth * hr.mort_mil))
    a.TP["mod_c",  "dth",] <- 1-exp(-(v.r.dth * hr.mort_mod))
    a.TP["sev_c",  "dth",] <- 1-exp(-(v.r.dth * hr.mort_sev))
    a.TP["dth"  ,  "dth",] <- 1
    a.TP["mci_i"  ,"dth",] <- 1-exp(-(v.r.dth * hr.mort_mci))
    a.TP["mil_i"  ,"dth",] <- 1-exp(-(v.r.dth * hr.mort_mil))
    a.TP["mod_i"  ,"dth",] <- 1-exp(-(v.r.dth * hr.mort_mod))
    a.TP["sev_i"  ,"dth",] <- 1-exp(-(v.r.dth * hr.mort_sev))
    
    # TP matrix state: from mci-on community-setting
    a.TP["mcion_c","mcion_c",] <- v.p.mcion_mci * (1-p.mci_i) * (1-v.p.discontinuation) * (1-a.TP["mcion_c","dth",])
    a.TP["mcion_c","mciof_c",] <- v.p.mcion_mci * (1-p.mci_i) *    v.p.discontinuation  * (1-a.TP["mcion_c","dth",])
    a.TP["mcion_c","milon_c",] <- v.p.mcion_mil * (1-p.mci_i) * (1-v.p.discontinuation) * (1-a.TP["mcion_c","dth",])
    a.TP["mcion_c","milof_c",] <- v.p.mcion_mil * (1-p.mci_i) *    v.p.discontinuation  * (1-a.TP["mcion_c","dth",])
    a.TP["mcion_c","mod_c",]   <- v.p.mcion_mod * (1-p.mci_i)                           * (1-a.TP["mcion_c","dth",])
    a.TP["mcion_c","sev_c",]   <- v.p.mcion_sev * (1-p.mci_i)                           * (1-a.TP["mcion_c","dth",])
    a.TP["mcion_c","mci_i",]   <- v.p.mcion_mci *    p.mci_i                            * (1-a.TP["mcion_c","dth",])
    a.TP["mcion_c","mil_i",]   <- v.p.mcion_mil *    p.mci_i                            * (1-a.TP["mcion_c","dth",])
    a.TP["mcion_c","mod_i",]   <- v.p.mcion_mod *    p.mci_i                            * (1-a.TP["mcion_c","dth",])
    a.TP["mcion_c","sev_i",]   <- v.p.mcion_sev *    p.mci_i                            * (1-a.TP["mcion_c","dth",])
    
    # TP matrix state: from mci-off community-setting
    a.TP["mciof_c","mciof_c",] <- v.p.mci_mci   * (1-p.mci_i)                           * (1-a.TP["mciof_c","dth",])
    a.TP["mciof_c","milof_c",] <- v.p.mci_mil   * (1-p.mci_i)                           * (1-a.TP["mciof_c","dth",])
    a.TP["mciof_c","mod_c",]   <- v.p.mci_mod   * (1-p.mci_i)                           * (1-a.TP["mciof_c","dth",])
    a.TP["mciof_c","sev_c",]   <- v.p.mci_sev   * (1-p.mci_i)                           * (1-a.TP["mciof_c","dth",])
    a.TP["mciof_c","mci_i",]   <- v.p.mci_mci   *    p.mci_i                            * (1-a.TP["mciof_c","dth",])
    a.TP["mciof_c","mil_i",]   <- v.p.mci_mil   *    p.mci_i                            * (1-a.TP["mciof_c","dth",])
    a.TP["mciof_c","mod_i",]   <- v.p.mci_mod   *    p.mci_i                            * (1-a.TP["mciof_c","dth",])
    a.TP["mciof_c","sev_i",]   <- v.p.mci_sev   *    p.mci_i                            * (1-a.TP["mciof_c","dth",])
    
    # TP matrix state: from mild-on community-setting
    a.TP["milon_c","mcion_c",] <- v.p.milon_mci * (1-p.mil_i) * (1-v.p.discontinuation) * (1-a.TP["milon_c","dth",])
    a.TP["milon_c","mciof_c",] <- v.p.milon_mci * (1-p.mil_i) *    v.p.discontinuation  * (1-a.TP["milon_c","dth",])
    a.TP["milon_c","milon_c",] <- v.p.milon_mil * (1-p.mil_i) * (1-v.p.discontinuation) * (1-a.TP["milon_c","dth",])
    a.TP["milon_c","milof_c",] <- v.p.milon_mil * (1-p.mil_i) *    v.p.discontinuation  * (1-a.TP["milon_c","dth",])
    a.TP["milon_c","mod_c",]   <- v.p.milon_mod * (1-p.mil_i)                           * (1-a.TP["milon_c","dth",])
    a.TP["milon_c","sev_c",]   <- v.p.milon_sev * (1-p.mil_i)                           * (1-a.TP["milon_c","dth",])
    a.TP["milon_c","mci_i",]   <- v.p.milon_mci *    p.mil_i                            * (1-a.TP["milon_c","dth",])
    a.TP["milon_c","mil_i",]   <- v.p.milon_mil *    p.mil_i                            * (1-a.TP["milon_c","dth",])
    a.TP["milon_c","mod_i",]   <- v.p.milon_mod *    p.mil_i                            * (1-a.TP["milon_c","dth",])
    a.TP["milon_c","sev_i",]   <- v.p.milon_sev *    p.mil_i                            * (1-a.TP["milon_c","dth",])
    
    # TP matrix state: from mild-off community-setting
    a.TP["milof_c","mciof_c",] <- v.p.mil_mci   * (1-p.mil_i)                           * (1-a.TP["milof_c","dth",])
    a.TP["milof_c","milof_c",] <- v.p.mil_mil   * (1-p.mil_i)                           * (1-a.TP["milof_c","dth",])
    a.TP["milof_c","mod_c",]   <- v.p.mil_mod   * (1-p.mil_i)                           * (1-a.TP["milof_c","dth",])
    a.TP["milof_c","sev_c",]   <- v.p.mil_sev   * (1-p.mil_i)                           * (1-a.TP["milof_c","dth",])
    a.TP["milof_c","mci_i",]   <- v.p.mil_mci   *    p.mil_i                            * (1-a.TP["milof_c","dth",])
    a.TP["milof_c","mil_i",]   <- v.p.mil_mil   *    p.mil_i                            * (1-a.TP["milof_c","dth",])
    a.TP["milof_c","mod_i",]   <- v.p.mil_mod   *    p.mil_i                            * (1-a.TP["milof_c","dth",])
    a.TP["milof_c","sev_i",]   <- v.p.mil_sev   *    p.mil_i                            * (1-a.TP["milof_c","dth",])
    
    # TP matrix state: from moderate community-setting
    a.TP["mod_c","milof_c",] <- v.p.mod_mil     * (1-p.mod_i)                           * (1-a.TP["mod_c","dth",])
    a.TP["mod_c","mod_c",]   <- v.p.mod_mod     * (1-p.mod_i)                           * (1-a.TP["mod_c","dth",])
    a.TP["mod_c","sev_c",]   <- v.p.mod_sev     * (1-p.mod_i)                           * (1-a.TP["mod_c","dth",])
    a.TP["mod_c","mil_i",]   <- v.p.mod_mil     *    p.mod_i                            * (1-a.TP["mod_c","dth",])
    a.TP["mod_c","mod_i",]   <- v.p.mod_mod     *    p.mod_i                            * (1-a.TP["mod_c","dth",])
    a.TP["mod_c","sev_i",]   <- v.p.mod_sev     *    p.mod_i                            * (1-a.TP["mod_c","dth",])
    
    # TP matrix state: from severe community-setting
    a.TP["sev_c","milof_c",] <- v.p.sev_mil     * (1-p.sev_i)                           * (1-a.TP["sev_c","dth",])
    a.TP["sev_c","mod_c",]   <- v.p.sev_mod     * (1-p.sev_i)                           * (1-a.TP["sev_c","dth",])
    a.TP["sev_c","sev_c",]   <- v.p.sev_sev     * (1-p.sev_i)                           * (1-a.TP["sev_c","dth",])
    a.TP["sev_c","mil_i",]   <- v.p.sev_mil     *    p.sev_i                            * (1-a.TP["sev_c","dth",])
    a.TP["sev_c","mod_i",]   <- v.p.sev_mod     *    p.sev_i                            * (1-a.TP["sev_c","dth",])
    a.TP["sev_c","sev_i",]   <- v.p.sev_sev     *    p.sev_i                            * (1-a.TP["sev_c","dth",])
    
    # TP matrix state: from mci institutionalized-setting
    a.TP["mci_i","mci_i",] <- v.p.mci_mci                                               * (1-a.TP["mci_i","dth",])
    a.TP["mci_i","mil_i",] <- v.p.mci_mil                                               * (1-a.TP["mci_i","dth",])
    a.TP["mci_i","mod_i",] <- v.p.mci_mod                                               * (1-a.TP["mci_i","dth",])
    a.TP["mci_i","sev_i",] <- v.p.mci_sev                                               * (1-a.TP["mci_i","dth",])
    
    # TP matrix state: from mild institutionalized-setting
    a.TP["mil_i","mci_i",] <- v.p.mil_mci                                               * (1-a.TP["mil_i","dth",])
    a.TP["mil_i","mil_i",] <- v.p.mil_mil                                               * (1-a.TP["mil_i","dth",])
    a.TP["mil_i","mod_i",] <- v.p.mil_mod                                               * (1-a.TP["mil_i","dth",])
    a.TP["mil_i","sev_i",] <- v.p.mil_sev                                               * (1-a.TP["mil_i","dth",])
    
    # TP matrix state: from moderate institutionalized-setting
    a.TP["mod_i","mil_i",] <- v.p.mod_mil                                               * (1-a.TP["mod_i","dth",])
    a.TP["mod_i","mod_i",] <- v.p.mod_mod                                               * (1-a.TP["mod_i","dth",])
    a.TP["mod_i","sev_i",] <- v.p.mod_sev                                               * (1-a.TP["mod_i","dth",])
    
    # TP matrix state: from severe institutionalized-setting
    a.TP["sev_i","mil_i",] <- v.p.sev_mil                                               * (1-a.TP["sev_i","dth",])
    a.TP["sev_i","mod_i",] <- v.p.sev_mod                                               * (1-a.TP["sev_i","dth",])
    a.TP["sev_i","sev_i",] <- v.p.sev_sev                                               * (1-a.TP["sev_i","dth",])  
    
    # check TPs are within 0-1 range
    if(any(a.TP<0 | a.TP>1)) stop("one or more transition probabilities are lower than 0 or higher than 1")
    # check TPs sum to 1 for each cycle (STEP B5: some checks)
    for(i in v.names_state) {
      temp1 <- colSums(a.TP[i,,])
      if(!isTRUE(all.equal(current = temp1, target = rep(1,n.cycle), tolerance = 1e-10))) stop(paste("TPs for",i,"do not add up to 1"))
    }
    
    # initialize state trace (STEP B6: initialize objects to store strategy outcomes)
    m.trace <- matrix(data = NA, nrow = n.cycle, ncol = n.state, dimnames = list(NULL,v.names_state))
    
    # initialize output table
    m.out <- matrix(data = NA, nrow = n.cycle, ncol = 29, dimnames = list(NULL, c(colnames(m.trace),"qaly_pt","qaly_ic","cost_dx","cost_tx","cost_hc","cost_sc","cost_ic","mci","mil","mod","sev","commun","instit","ontx","ly","qaly","cost","nhb")))
    
    # set starting state distribution (STEP B7: starting state)
    m.trace[1,] <- m.trace1
    
    # markov multiplication (STEP B8: markov multiplication by looping over cycles)
    for(t in 1:(n.cycle-1)) {
      m.trace[t+1,] <- m.trace[t,] %*% a.TP[,,t]
    }
    
    # primary economic outputs (STEP B9: multiply states with utility and cost estimates)
    m.out[,colnames(m.trace)] <- m.trace
    m.out[,"qaly_pt"] <- m.trace %*% c(u.mci_pt, u.mci_pt, u.mil_pt, u.mil_pt, u.mod_pt, u.sev_pt, u.mci_pt_i, u.mil_pt_i, u.mod_pt_i, u.sev_pt_i, 0) # must match order of states
    m.out[,"qaly_ic"] <- m.trace %*% c(u.mci_ic, u.mci_ic, u.mil_ic, u.mil_ic, u.mod_ic, u.sev_ic, u.mci_ic_i, u.mil_ic_i, u.mod_ic_i, u.sev_ic_i, 0) # must match order of states
    m.out[,"cost_dx"] <- 0
    m.out[,"cost_tx"] <- m.trace %*% c(c.Tx    , 0       , c.Tx    , 0       , 0       , 0       , 0         , 0         , 0         , 0         , 0) # must match order of states
    m.out[,"cost_hc"] <- m.trace %*% c(c.mci_hc, c.mci_hc, c.mil_hc, c.mil_hc, c.mod_hc, c.sev_hc, c.mci_hc_i, c.mil_hc_i, c.mod_hc_i, c.sev_hc_i, 0) # must match order of states
    m.out[,"cost_sc"] <- m.trace %*% c(c.mci_sc, c.mci_sc, c.mil_sc, c.mil_sc, c.mod_sc, c.sev_sc, c.mci_sc_i, c.mil_sc_i, c.mod_sc_i, c.sev_sc_i, 0) # must match order of states
    m.out[,"cost_ic"] <- m.trace %*% c(c.mci_ic, c.mci_ic, c.mil_ic, c.mil_ic, c.mod_ic, c.sev_ic, c.mci_ic_i, c.mil_ic_i, c.mod_ic_i, c.sev_ic_i, 0) # must match order of states
    
    # half-cycle correction (STEP B10: apply half-cycle correction to all outcomes)
    if(half_cycle_correction) {
      for (j in colnames(m.out)) {
        for (i in 1:(n.cycle-1)) {
          m.out[i,j]   <- (m.out[i,j]   + m.out[i+1,j])   * 0.5
        }
      }
      m.out <- m.out[-n.cycle,] # remove the last cycle
    }
    
    # add additional inputs to cycle 1
    if(strat=="int") {
      m.out[,"qaly_pt"][1] <- m.out[,"qaly_pt"][1] + u.Tx_start
      m.out[,"cost_dx"][1] <- m.out[,"cost_dx"][1] + c.Tx_start
    }
    
    # define vector for discounting QALYs, costs and effects (STEP B11: apply discounting to all outcomes)
    n <- ifelse(test=half_cycle_correction, yes=2, no=1)
    v.discount_EFFECT <- 1 / (( 1 + discount_EFFECT) ^ (0 : (n.cycle-n)))
    v.discount_QALY   <- 1 / (( 1 + discount_QALY)   ^ (0 : (n.cycle-n)))
    v.discount_COST   <- 1 / (( 1 + discount_COST)   ^ (0 : (n.cycle-n)))
    
    # apply discounting to all outcomes
    for(i in c(colnames(m.trace),"ly")) {
      m.out[,i] <- m.out[,i]*v.discount_EFFECT
    }
    for(i in c("qaly_pt","qaly_ic")) {
      m.out[,i] <- m.out[,i]*v.discount_QALY
    }
    for(i in c("cost_dx","cost_tx","cost_hc","cost_sc","cost_ic")) {
      m.out[,i] <- m.out[,i]*v.discount_COST
    }
    
    # totals
    m.out[,"mci"] <- rowSums(m.out[,c("mcion_c","mciof_c","mci_i")])
    m.out[,"mil"] <- rowSums(m.out[,c("milon_c","milof_c","mil_i")])
    m.out[,"mod"] <- rowSums(m.out[,c("mod_c","mod_i")])
    m.out[,"sev"] <- rowSums(m.out[,c("sev_c","sev_i")])
    m.out[,"commun"] <- rowSums(m.out[,c("mcion_c","mciof_c","milon_c","milof_c","mod_c","sev_c")])
    m.out[,"instit"] <- rowSums(m.out[,c("mci_i","mil_i","mod_i","sev_i")])
    m.out[,"ontx"] <- rowSums(m.out[,c("mcion_c","milon_c")])
    m.out[,"qaly"] <- rowSums(m.out[,c("qaly_pt","qaly_ic")])
    m.out[,"cost"] <- rowSums(m.out[,c("cost_dx","cost_tx","cost_hc","cost_sc","cost_ic")])
    m.out[,"ly"] <- rowSums(m.out[,c("mci","mil","mod","sev")])
    
    # calculate net health benefit
    m.out[,"nhb"] <- m.out[,"qaly"] - (m.out[,"cost"] / wtp)
    
    # store strategy-specific output (STEP B12: store outcomes to be wrapped up by the 'run scenario' function)
    return(list(
      a.TP = a.TP, 
      m.trace = m.trace, 
      m.out = m.out
    ))
  }
  )
}



######################################## 2.2. RUN SCENARIO (STEP A1-A5) ########################################

# run scenario (STEP A: function to run a a scenario, which includes a loop over all strategies)
f.run_scenario <- function(l.inputs, detailed=FALSE) {

  # store strategy names and count (STEP A1: initialize objects to store scenario and strategy outcomes)
  v.names_strat <- l.inputs[["v.names_strat"]]
  n.strat <- length(v.names_strat) # number of strategies
  
  # initialize data frame to store scenario outcomes
  df.out <- data.frame(
    strategy = v.names_strat,
    QALY = numeric(n.strat),
    COST = numeric(n.strat),
    LY = numeric(n.strat),
    NHB = numeric(n.strat),
    row.names = v.names_strat, 
    stringsAsFactors = FALSE
  )
  
  # initialize list to store strategy outcomes (create empty list to store outcomes of each strategy)
  l.out <- vector(mode = "list", length = 0)
  
  # loop over strategies (STEP A2: run each strategy in a loop)
  for(strat in v.names_strat) {
    
    # run strategy (STEP B: run the strategy)
    l.strat <- f.run_strategy(l.inputs = l.inputs, strat = strat)
    
    # store strategy-specific output (STEP A3 store strategy results)
    l.out[[strat]] <- l.strat
    
    # store output (STEP A4 add strategy results to scenario outcomes)
    m.out <- l.strat[["m.out"]] # extract cycle-specific outcomes
    df.out[strat,"strategy"] <- strat # store strategy name
    df.out[strat,"QALY"] <- sum(m.out[,"qaly"]) # calculate total QALYs and store them
    df.out[strat,"COST"] <- sum(m.out[,"cost"]) # calculate total costs and store them
    df.out[strat,"LY"]   <- sum(m.out[,"ly"]) # calculate total live years and store them
    df.out[strat,"NHB"]  <- sum(m.out[,"nhb"]) # calculate total NHB and store them
  }
  
  # return basic result
  if(!detailed) return(df.out)
  
  # return detailed result
  if(detailed) {
    return(list(
      df.out = df.out,
      l.out = l.out
    ))
  }
}



######################################## 3. MODEL CALIBRATION ########################################

# n/a



######################################## 4. VALIDATION ########################################

# n/a



######################################## 5. ANALYSIS ########################################

######################################## 5.1. CROSS-VALIDATION: ICER ########################################

if(F) {
  
  # run scenario
  l.out_icer <- f.run_scenario(l.inputs = l.inputs_icer, detailed = TRUE)
  
  # additional results
  m.result_icer <- matrix(data = NA, nrow = ncol(l.out_icer$l.out$soc$m.out), ncol = 3, dimnames = list(colnames(l.out_icer$l.out$soc$m.out), c("soc","int","dif")))
  m.result_icer[,"soc"] <- colSums(l.out_icer$l.out$soc$m.out)
  m.result_icer[,"int"] <- colSums(l.out_icer$l.out$int$m.out)
  m.result_icer[,"dif"] <- m.result_icer[,"int"] - m.result_icer[,"soc"]
  
  # compare to publication
  print(round(m.result_icer[c("ly","qaly","cost"),c("soc","int","dif")],2))
  icer_icer <- calculate_icers(cost = l.out_icer[["df.out"]][,"COST"], effect = l.out_icer[["df.out"]][,"QALY"], strategies = l.out_icer[["df.out"]][,"strategy"])
  print(icer_icer)
  
  # proportion in state
  l.inputs_icer_nohccdis <- l.inputs_icer # run model without half-cycle correction and discounted effects
  l.inputs_icer_nohccdis[["discount_EFFECT"]] <- 0
  l.inputs_icer_nohccdis[["half_cycle_correction"]] <- FALSE
  l.out_icer_nohccdis <- f.run_scenario(l.inputs = l.inputs_icer_nohccdis, detailed = TRUE)
    # write.table(x = round(l.out_icer_nohccdis$l.out$soc$m.out[1:11,c("mci","mil","mod","sev","dth")],2), file = "clipboard", sep = "\t", row.names = FALSE, col.names = FALSE)
    # write.table(x = round(l.out_icer_nohccdis$l.out$int$m.out[1:11,c("mci","mil","mod","sev","dth")],2), file = "clipboard", sep = "\t", row.names = FALSE, col.names = FALSE)
  
  # export results for IPECAD repository
    # write.table(x = m.result_icer[c("cost_hc","cost_sc","cost_ic","cost_tx","qaly_pt","qaly_ic"),c("soc","int")], file = "clipboard", sep = "\t", row.names = FALSE, col.names = FALSE)
  l.inputs_icer_repository <- l.inputs_icer # run model without half-cycle correction and discounted effects
  l.inputs_icer_repository[["half_cycle_correction"]] <- FALSE
  l.inputs_icer_repository[["discount_EFFECT"]] <- 0
  l.out_icer_repository <- f.run_scenario(l.inputs = l.inputs_icer_repository, detailed = TRUE)
    # write.table(x = l.out_icer_repository$l.out$soc$m.out[1:26,c("mci","mil","mod","sev","dth")], file = "clipboard", sep = "\t", row.names = FALSE, col.names = FALSE)
    # write.table(x = l.out_icer_repository$l.out$int$m.out[1:26,c("mci","mil","mod","sev","dth")], file = "clipboard", sep = "\t", row.names = FALSE, col.names = FALSE)
  
  # standard tables/plots
  
  if(T) {
    
    # all outcomes
    # str(l.out_icer) # show structure of output
    # l.out_icer[["df.out"]]
    # l.out_icer # all detailed output
    # l.out_icer[["l.out"]] # strategy details
    # l.out_icer[["l.out"]][["soc"]] # strategy 'soc' details
    # l.out_icer[["l.out"]][["int"]] # strategy 'int' details
    # l.out_icer[["l.out"]][["soc"]][["a.TP"]] # strategy 'soc' transition probability matrix
    # l.out_icer[["l.out"]][["soc"]][["m.trace"]] # strategy 'soc' state trace
    # l.out_icer[["l.out"]][["soc"]][["m.out"]] # strategy 'soc' outcomes per cycle
    # l.out_icer[["l.out"]][["int"]][["a.TP"]] # idem for 'int'
    # l.out_icer[["l.out"]][["int"]][["m.trace"]]
    # l.out_icer[["l.out"]][["int"]][["m.out"]]
    
    # table: summary
    table.summary_data <- m.result_icer[c("mci","mil","mod","sev","ly","qaly","cost"),c("soc","int","dif")]
    table.summary <- format(
      table.summary_data, 
      digits=2, 
      scientific=FALSE, 
      big.mark=","
    )
    
    print(round(table.summary_data,2))
    print(table.summary)
    
    # plot: time in state
    plot.timestate_data <- m.result_icer[c("mci","mil","mod","sev"),c("int","soc")]
    
      # windows(width=7, height=4, pointsize=12)
    par(mar=c(5, 4, 4, 1), xpd=TRUE)
    barplot(
      height = plot.timestate_data, 
      horiz = TRUE, 
      xaxt = "n", 
      xlim = c(0,ceiling(max(colSums(plot.timestate_data)))), 
      xlab = "time (years)", 
      ylab = "strategy", 
      col=c("green","yellow","orange","red"), 
      space = 0.2, 
      main = "mean time in state"
    )
    axis(side = 1, at = 0:ceiling(max(colSums(plot.timestate_data))))
    text(x=c(0,cumsum(plot.timestate_data[1:3,"int"])), y=1, labels=round(plot.timestate_data[,"int"],1), pos=4)
    text(x=c(0,cumsum(plot.timestate_data[1:3,"soc"])), y=2, labels=round(plot.timestate_data[,"soc"],1), pos=4)
    text(x=c(0,cumsum(plot.timestate_data[1:3,"int"])), y=0.5, labels=c("MCI","mild","moderate","severe"), pos=4)
    #legend(x="bottom", legend=c("mci","mil","mod","sev"), inset=c(0,-0.6), horiz=TRUE, fill=c("green","yellow","orange","red"))
    
    # table: state trace
    tableplot.statetracesoc_data <- cbind(
      mci = rowSums(l.out_icer[["l.out"]][["soc"]][["m.trace"]][,c("mcion_c","mciof_c","mci_i")]), 
      mil = rowSums(l.out_icer[["l.out"]][["soc"]][["m.trace"]][,c("milon_c","milof_c","mil_i")]),
      mod = rowSums(l.out_icer[["l.out"]][["soc"]][["m.trace"]][,c("mod_c","mod_i")]), 
      sev = rowSums(l.out_icer[["l.out"]][["soc"]][["m.trace"]][,c("sev_c","sev_i")]), 
      dth =         l.out_icer[["l.out"]][["soc"]][["m.trace"]][,c("dth")]
    )
    tableplot.statetraceint_data <- cbind(
      mci = rowSums(l.out_icer[["l.out"]][["int"]][["m.trace"]][,c("mcion_c","mciof_c","mci_i")]), 
      mil = rowSums(l.out_icer[["l.out"]][["int"]][["m.trace"]][,c("milon_c","milof_c","mil_i")]),
      mod = rowSums(l.out_icer[["l.out"]][["int"]][["m.trace"]][,c("mod_c","mod_i")]), 
      sev = rowSums(l.out_icer[["l.out"]][["int"]][["m.trace"]][,c("sev_c","sev_i")]), 
      dth =         l.out_icer[["l.out"]][["int"]][["m.trace"]][,c("dth")]
    )
    
    print(round(tableplot.statetracesoc_data[1:min(nrow(tableplot.statetracesoc_data),10),],2))
    print(round(tableplot.statetraceint_data[1:min(nrow(tableplot.statetracesoc_data),10),],2))
    
    # plot: state trace
    v.age_range <- c(l.inputs_icer[["age_start"]]:(l.inputs_icer[["age_start"]]+l.inputs_icer[["n.cycle"]]-1)) # store age range
    xx <- c(v.age_range, rev(v.age_range)) # prepare polygon x-values
    yy_mci <- c(tableplot.statetracesoc_data[,"mci"], rev(tableplot.statetraceint_data[,"mci"])) # polygon y-values
    yy_mil <- c(tableplot.statetracesoc_data[,"mil"], rev(tableplot.statetraceint_data[,"mil"])) # idem
    yy_mod <- c(tableplot.statetracesoc_data[,"mod"], rev(tableplot.statetraceint_data[,"mod"])) # idem
    yy_sev <- c(tableplot.statetracesoc_data[,"sev"], rev(tableplot.statetraceint_data[,"sev"])) # idem
    yy_dth <- c(tableplot.statetracesoc_data[,"dth"], rev(tableplot.statetraceint_data[,"dth"])) # idem
    
      # windows(width=7, height=7, pointsize=12)
    par(mar=c(5, 4, 4, 1)+0.1, xpd=FALSE)
    matplot(
      x = v.age_range, 
      y = cbind(tableplot.statetracesoc_data,tableplot.statetraceint_data), 
      type = "n", 
      xlab = "age", 
      ylab = "proportion in state", 
      ylim = c(0,1), 
      main = "state trace"
    )
    polygon(xx, yy_mci, col = "gray95", border = FALSE)
    polygon(xx, yy_mil, col = "gray95", border = FALSE)
    polygon(xx, yy_mod, col = "gray95", border = FALSE)
    polygon(xx, yy_sev, col = "gray95", border = FALSE)
    polygon(xx, yy_dth, col = "gray95", border = FALSE)
    matlines(
      x = v.age_range, 
      y = tableplot.statetracesoc_data, 
      type = "l",
      lty = 1,
      col = c("green","yellow","orange","red","black")
    )
    matlines(
      x = v.age_range, 
      y = tableplot.statetraceint_data, 
      type = "l",
      lty = 2,
      col = c("green","yellow","orange","red","black")
    )
    legend(x="topright", legend=c("mci","mil","mod","sev","dth"), col=c("green","yellow","orange","red","black"), lty=1, bg="white")
    legend(x="right", legend=c("soc","int"), col="black", lty=c(1,2), bg="white")
    
    # table: incremental cost-effectiveness ratio
    print(as.data.frame(icer_icer))
    
    # plot: incremental cost-effectiveness plane
    par(mar=c(5, 4, 4, 1)+0.1, xpd=FALSE)
      # windows(width=7, height=7, pointsize=12)
    print(plot(icer_icer, label="all"))
    
    # plot: cost difference by sector over time
    m.cost_incr_pos <- m.cost_incr_neg <- l.out_icer[["l.out"]][["int"]][["m.out"]][,c("cost_dx","cost_tx","cost_hc","cost_sc","cost_ic")] - l.out_icer[["l.out"]][["soc"]][["m.out"]][,c("cost_dx","cost_tx","cost_hc","cost_sc","cost_ic")] # split positive and negative
    m.cost_incr_pos[m.cost_incr_pos<0] <- 0
    m.cost_incr_neg[m.cost_incr_neg>=0] <- 0
    
      # windows(width=7, height=7, pointsize=12)
    par(mar=c(5, 4, 4, 2)+0.1, xpd=FALSE)
    barplot(
      height = t(m.cost_incr_pos),
      beside = F,
      xlab = "time (years)",
      ylab = "annual incremental costs",
      ylim = c(
        min(m.cost_incr_neg) + min(m.cost_incr_neg)*0.10, 
        max(m.cost_incr_pos) + max(m.cost_incr_pos)*0.10
      ),
      col = rainbow(5), 
      names.arg = 1:nrow(m.cost_incr_pos),
      main = "costs by sector over time"
    )
    barplot(
      height = t(m.cost_incr_neg),
      beside = F,
      col = rainbow(5),
      add = T
    )
    legend(x = "topright", legend = c("diagnostic","treatment","health","social","informal"), fill = rainbow(5))
    
  }
  
}



######################################## 5.2. CROSS-VALIDATION: AD-ACE ########################################

if(F) {
  
  # run scenario and results
  l.out_adace <- f.run_scenario(l.inputs = l.inputs_adace, detailed = TRUE)
  
  # additional results
  m.result_adace <- matrix(data = NA, nrow = ncol(l.out_adace$l.out$soc$m.out), ncol = 3, dimnames = list(colnames(l.out_adace$l.out$soc$m.out), c("soc","int","dif")))
  m.result_adace[,"soc"] <- colSums(l.out_adace$l.out$soc$m.out)
  m.result_adace[,"int"] <- colSums(l.out_adace$l.out$int$m.out)
  m.result_adace[,"dif"] <- m.result_adace[,"int"] - m.result_adace[,"soc"]
  
  # compare to publication
  print(round(m.result_adace[c("ly","qaly","cost"),c("soc","int","dif")],2))
  icer_adace <- calculate_icers(cost = l.out_adace[["df.out"]][,"COST"], effect = l.out_adace[["df.out"]][,"QALY"], strategies = l.out_adace[["df.out"]][,"strategy"])
  print(icer_adace)
  
  # proportion in state
  l.inputs_adace_nohccdis <- l.inputs_adace # run model without half-cycle correction and discounted effects
  l.inputs_adace_nohccdis[["discount_EFFECT"]] <- 0
  l.inputs_adace_nohccdis[["half_cycle_correction"]] <- FALSE
  l.out_adace_nohccdis <- f.run_scenario(l.inputs = l.inputs_adace_nohccdis, detailed = TRUE)
    # write.table(x = round(l.out_adace_nohccdis$l.out$soc$m.out[1:11,c("mci","mil","mod","sev","dth")],2), file = "clipboard", sep = "\t", row.names = FALSE, col.names = FALSE)
    # write.table(x = round(l.out_adace_nohccdis$l.out$int$m.out[1:11,c("mci","mil","mod","sev","dth")],2), file = "clipboard", sep = "\t", row.names = FALSE, col.names = FALSE)
  
  # export results for IPECAD repository
    # write.table(x = m.result_adace[c("cost_hc","cost_sc","cost_ic","cost_tx","qaly_pt","qaly_ic"),c("soc","int")], file = "clipboard", sep = "\t", row.names = FALSE, col.names = FALSE)
  l.inputs_adace_repository <- l.inputs_adace # run model without half-cycle correction and discounted effects
  l.inputs_adace_repository[["half_cycle_correction"]] <- FALSE
  l.inputs_adace_repository[["discount_EFFECT"]] <- 0
  l.out_adace_repository <- f.run_scenario(l.inputs = l.inputs_adace_repository, detailed = TRUE)
    # write.table(x = l.out_adace_repository$l.out$soc$m.out[1:26,c("mci","mil","mod","sev","dth")], file = "clipboard", sep = "\t", row.names = FALSE, col.names = FALSE)
    # write.table(x = l.out_adace_repository$l.out$int$m.out[1:26,c("mci","mil","mod","sev","dth")], file = "clipboard", sep = "\t", row.names = FALSE, col.names = FALSE)
  
}



######################################## 5.3. UNCERTAINTY SCENARIOS ########################################

if(F) {
  
  # base case
  # icer (already defined)
  
  # CDR-SB RR=23%
  l.inputs_icer_2 <- l.inputs_icer
  l.inputs_icer_2[["rr.tx_mci_mil"]] <- 1-0.23
  l.inputs_icer_2[["rr.tx_mil_mod"]] <- 1-0.23
  l.inputs_icer_2[["rr.tx_mil_sev"]] <- 1-0.23
  
  # CDR-SB time shift (calibrate)
  # see code chapter 'calibrate to time shift'
  
  # MMSE progression (Vos & SveDem)
  ## source: [Vos, 2015: https://doi.org/10.1093/brain/awv029]
  ## option: Amyloid positive & neuronal loss undetermined
  ### Operationalized by diagnostic criteria NIA-AA categories: 'NIA-AA high AD' (Amyloid+, Injury+) and 'conflicting IAP' (Amyloid+, Injury-)
  ### corresponding 3-year cumulative incidence probability: 'high AD' = 59% (AD dementia; table 3) and 4% (non-AD dementia; mentioned in text), and 22% (AD dementia; table 3) and 4% (non-AD dementia; mentioned in text) with prevalence of 353 and 49 respectively (respectively)
  ### This results into a weighted 3-year cumulative incidence of (converting all 4 probabilities to rates before weighting and averaging):
  temp.est2 <- 1-exp(- ( (-log(1-0.59) + -log(1-0.04))*353 + (-log(1-0.22) + -log(1-0.04))*49 ) / (353+49) )
  ## and corresponding 1-year probability of 
  temp.est2 <- 1-exp(- -log(1-temp.est2)/3)
  temp.est2
  ## alternative option: Amyloid positive & neuronal loss positive
  ### Operationalized by diagnostic criteria NIA-AA categories: 'NIA-AA high AD' (Amyloid+, Injury+)
  ### corresponding 3-year cumulative incidence probability: 'high AD' = 59% (AD dementia; table 3) and 4% (non-AD dementia; mentioned in text)
  ### This results into a weighted 3-year cumulative incidence of:
  temp.est1 <- 1-exp(- (-log(1-0.59) + -log(1-0.04)) )
  ## and corresponding 1-year probability of 
  temp.est1 <- 1-exp(- -log(1-temp.est1)/3)
  temp.est1
  
  l.inputs_icer_4 <- l.inputs_icer
  l.inputs_icer_4[["p.mci_mil"]] <- round(temp.est2,3) # 0.248
  l.inputs_icer_4[["p.mci_mod"]] <- 0
  l.inputs_icer_4[["p.mci_sev"]] <- 0
  l.inputs_icer_4[["p.mil_mci"]] <- 0
  l.inputs_icer_4[["p.mil_mod"]] <- 0.293
  l.inputs_icer_4[["p.mil_sev"]] <- 0.001
  l.inputs_icer_4[["p.mod_mil"]] <- 0.087
  l.inputs_icer_4[["p.mod_sev"]] <- 0.109
  l.inputs_icer_4[["p.sev_mil"]] <- 0.000
  l.inputs_icer_4[["p.sev_mod"]] <- 0.196
  
  # Mortality (SveDem)
  l.inputs_icer_5 <- l.inputs_icer
  l.inputs_icer_5[["hr.mort_mci"]] <- 1
  l.inputs_icer_5[["hr.mort_mil"]] <- 1.318 * 1.82
  l.inputs_icer_5[["hr.mort_mod"]] <- 2.419 * 1.82
  l.inputs_icer_5[["hr.mort_sev"]] <- 4.267 * 1.82
  
  # MMSE progression and mortality (Vos & SveDem)
  l.inputs_icer_6 <- l.inputs_icer_4
  l.inputs_icer_6[["hr.mort_mci"]] <- 1
  l.inputs_icer_6[["hr.mort_mil"]] <- 1.318 * 1.82
  l.inputs_icer_6[["hr.mort_mod"]] <- 2.419 * 1.82
  l.inputs_icer_6[["hr.mort_sev"]] <- 4.267 * 1.82
  
  # A: continue treatment & no waning during treatment
  # n/a: identical to base case
  
  # B: continue treatment & waning during treatment
  l.inputs_icer_B <- l.inputs_icer
  0.69^(1-0.30)^5 # 30% waning reduces RR to 0.94 (compared to 0.69) after 5 years
  l.inputs_icer_B[["tx_waning"]] <- 0.30
  
  # C: stop treatment & no waning after treatment stop
  l.inputs_icer_C <- l.inputs_icer
  l.inputs_icer_C[["rr.tx_mci_mil_dis"]] <- l.inputs_icer_C[["rr.tx_mci_mil"]]
  l.inputs_icer_C[["rr.tx_mil_mod_dis"]] <- l.inputs_icer_C[["rr.tx_mil_mod"]]
  l.inputs_icer_C[["rr.tx_mil_sev_dis"]] <- l.inputs_icer_C[["rr.tx_mil_sev"]]
  l.inputs_icer_C[["p.tx_discontinuation2"]] <- 1
  l.inputs_icer_C[["tx_discontinuation2_begin"]] <- 2
  
  # D: stop treatment & waning after treatment stop
  l.inputs_icer_D <- l.inputs_icer_C
  l.inputs_icer_D[["tx_waning_dis"]] <- 0.30
  
  # E: stop treatment & no effect after treatment stop
  l.inputs_icer_E <- l.inputs_icer
  l.inputs_icer_E[["p.tx_discontinuation2"]] <- 1
  l.inputs_icer_E[["tx_discontinuation2_begin"]] <- 2
  
  # methodology: cycle lengths 
  # see code chapter 'cycle time'
  
  # run scenarios
  
  ## initialize outcomes lists and tables
  l.out_scenario <- list()
  l.out_scenario_prep <- list()
  m.table1 <- matrix(data = NA, nrow = 12, ncol = 7, dimnames = list(NULL,c("mcimil","ly","qaly","cost_dxtx","cost_care","nhb","icer")))
  
  ## run scenarios
  l.out_scenario[[1]] <- f.run_scenario(l.inputs = l.inputs_icer, detailed = TRUE)
  l.out_scenario[[2]] <- f.run_scenario(l.inputs = l.inputs_icer_2, detailed = TRUE)
  
  l.out_scenario[[4]] <- f.run_scenario(l.inputs = l.inputs_icer_4, detailed = TRUE)
  l.out_scenario[[5]] <- f.run_scenario(l.inputs = l.inputs_icer_5, detailed = TRUE)
  l.out_scenario[[6]] <- f.run_scenario(l.inputs = l.inputs_icer_6, detailed = TRUE)
  
  l.out_scenario[[8]] <- f.run_scenario(l.inputs = l.inputs_icer_B, detailed = TRUE)
  l.out_scenario[[9]] <- f.run_scenario(l.inputs = l.inputs_icer_C, detailed = TRUE)
  l.out_scenario[[10]] <- f.run_scenario(l.inputs = l.inputs_icer_D, detailed = TRUE)
  l.out_scenario[[11]] <- f.run_scenario(l.inputs = l.inputs_icer_E, detailed = TRUE)
  
  ## store results
  for(i in c(1,2, 4,5,6, 8,9,10,11)) {
    # additional results
    l.out_scenario_prep[[i]] <- matrix(data = NA, nrow = ncol(l.out_scenario[[i]]$l.out$soc$m.out), ncol = 3, dimnames = list(colnames(l.out_scenario[[i]]$l.out$soc$m.out), c("soc","int","dif")))
    l.out_scenario_prep[[i]][,"soc"] <- colSums(l.out_scenario[[i]]$l.out$soc$m.out)
    l.out_scenario_prep[[i]][,"int"] <- colSums(l.out_scenario[[i]]$l.out$int$m.out)
    l.out_scenario_prep[[i]][,"dif"] <- l.out_scenario_prep[[i]][,"int"] - l.out_scenario_prep[[i]][,"soc"]
    # checks
    if(l.out_scenario_prep[[i]]["mci","dif"]<0) stop("lower person years")
    if(l.out_scenario_prep[[i]]["mil","dif"]<0) stop("lower person years")
    # select from additional results
    m.table1[i,"mcimil"] <- sum(l.out_scenario_prep[[i]][c("mci","mil"),"dif"])
    m.table1[i,"ly"] <- l.out_scenario_prep[[i]]["ly","dif"]
    m.table1[i,"qaly"] <- l.out_scenario_prep[[i]]["qaly","dif"]
    m.table1[i,"cost_dxtx"] <- sum(l.out_scenario_prep[[i]][c("cost_dx","cost_tx"),"dif"])
    m.table1[i,"cost_care"] <- sum(l.out_scenario_prep[[i]][c("cost_hc","cost_sc","cost_ic"),"dif"])
    m.table1[i,"nhb"] <- sum(l.out_scenario_prep[[i]]["nhb","dif"])
    m.table1[i,"icer"] <- calculate_icers(cost = l.out_scenario[[i]][["df.out"]][,"COST"], effect = l.out_scenario[[i]][["df.out"]][,"QALY"], strategies = l.out_scenario[[i]][["df.out"]][,"strategy"])[2,"ICER"]
  }
  print(round(m.table1,2))
    # write.table(x = round(m.table1,2), file = "clipboard", sep = "\t", row.names = FALSE, col.names = FALSE)
  
  # proportion results extrapolated
  ## store cycle-specific outcomes
  m.out_icer_soc <- l.out_icer$l.out$soc$m.out
  m.out_icer_int <- l.out_icer$l.out$int$m.out
  ## initialize matrix for summary outcomes
  m.result_icer <- matrix(
    data = NA, 
    nrow = ncol(m.out_icer_soc), 
    ncol = 11, 
    dimnames = list(
      colnames(m.out_icer_soc),
      c("soc","int","dif","dif_relative","soc_within","int_within","soc_extrapolate","int_extrapolate","dif_within","dif_extrapolate","dif_p_extrapolate")
    )
  )
  ## define number of cycles considered within-trial
  n.within <- 2
  ## fill matrix
  for(i in rownames(m.result_icer)) {
    m.result_icer[i,"soc"] <- sum(m.out_icer_soc[,i])
    m.result_icer[i,"int"] <- sum(m.out_icer_int[,i])
    m.result_icer[i,"soc_within"] <- sum(m.out_icer_soc[1:n.within,i])
    m.result_icer[i,"int_within"] <- sum(m.out_icer_int[1:n.within,i])
    m.result_icer[i,"soc_extrapolate"] <- sum(m.out_icer_soc[(n.within+1):nrow(m.out_icer_soc),i])
    m.result_icer[i,"int_extrapolate"] <- sum(m.out_icer_int[(n.within+1):nrow(m.out_icer_soc),i])
  }
  ## add absolute, relative and proportional difference
  m.result_icer[,"dif"] <- m.result_icer[,"int"] - m.result_icer[,"soc"]
  m.result_icer[,"dif_relative"] <- m.result_icer[,"dif"] / m.result_icer[,"soc"]
  m.result_icer[,"dif_within"] <- m.result_icer[,"int_within"] - m.result_icer[,"soc_within"]
  m.result_icer[,"dif_extrapolate"] <- m.result_icer[,"int_extrapolate"] - m.result_icer[,"soc_extrapolate"]
  m.result_icer[,"dif_p_extrapolate"] <- m.result_icer[,"dif_extrapolate"] / m.result_icer[,"dif"]
  ## select results for manuscript
  round(m.result_icer[c("mci","mil","ly","qaly","cost_dx","cost_tx","cost_hc","cost_sc","cost_ic"),c("dif_within","dif_extrapolate","dif_p_extrapolate")], digits=2) # note: proportional difference is invalid for outcomes with negative values
  
}



######################################## 5.3.1. CYCLE TIME ########################################

if(F) {
  
  # cycle time adjustment following guidance by Gidwani et al. [2020: https://doi.org/10.1007/s40273-020-00937-z]
  
  # model inputs list
  l.inputs_cycle <- l.inputs_icer
  
  # function to convert probability (p) to a different (cycle) time period (t=time proportional to current 1-year cycle length)
  f.p_time <- function(p, t) 1-(1-p)^(t)
  
  # proportion of the 1-year cycle length
  t <- 1/24 # half-month cycle length
  
  # simultaneously convert matrix of dementia transition probabilities to new cycle length
  ## temporary matrix of dementia transition probabilities
  t.TP <- matrix(data=NA, nrow=3, ncol=3, dimnames=list(c("mil","mod","sev"),c("mil","mod","sev")))
  t.TP["mil","mod"] <- l.inputs_cycle[["p.mil_mod"]]
  t.TP["mil","sev"] <- l.inputs_cycle[["p.mil_sev"]]
  t.TP["mod","mil"] <- l.inputs_cycle[["p.mod_mil"]]
  t.TP["mod","sev"] <- l.inputs_cycle[["p.mod_sev"]]
  t.TP["sev","mil"] <- l.inputs_cycle[["p.sev_mil"]]
  t.TP["sev","mod"] <- l.inputs_cycle[["p.sev_mod"]]
  t.TP[1,1] <- 1-sum(t.TP[1,-1])
  t.TP[2,2] <- 1-sum(t.TP[2,-2])
  t.TP[3,3] <- 1-sum(t.TP[3,-3])
  t.TP
  
  ## convert TPs using eigenvalues
  D <- diag(x=eigen(t.TP)$values); D
  V <- eigen(t.TP)$vectors; V
  t.TP2 <- V %*% D^(t) %*% solve(V); t.TP2
  ## manually check for non-numbers or negative values (if very small, might be rounded to 0)
  t.TP2
  t.TP2[t.TP2 < 0] <- 0 # convert small negative values to 0
  t.TP2 <- round(t.TP2, 3) # round
  dimnames(t.TP2) <- dimnames(t.TP) # add names
  rowSums(t.TP2) # check if rowsums add up to 1
  diag(t.TP2) <- NA # remove matrix diagonal (probability of remaining in the same state)
  t.TP2[1,1] <- 1-sum(t.TP2[1,-1]) # calculate probability of remaining in the same state as '1-(other probabilities)'
  t.TP2[2,2] <- 1-sum(t.TP2[2,-2])
  t.TP2[3,3] <- 1-sum(t.TP2[3,-3])
  t.TP2
  
  # check: convert back to 1-year cycle to estimate the difference with the original estimates
  D2 <- diag(x=eigen(t.TP2)$values); D2
  V2 <- eigen(t.TP2)$vectors; V2
  t.TP22 <- V2 %*% D2^(1/t) %*% solve(V2); t.TP22
  round(t.TP22,2)
  round(t.TP,2)
  round(t.TP - t.TP22,2)
  
  # check: compare p.mci_mil handled separately versus handled simultaneously with dementia transition probabilities (backtransitions from dementia to mci are assumed 0)
  ## handled separately
  t.TP2_mci <- rbind(0,cbind(0,t.TP2))
  rownames(t.TP2_mci)[1] <- "mci"
  colnames(t.TP2_mci)[1] <- "mci"
  t.TP2_mci["mci","mil"] <- f.p_time(p = l.inputs_cycle[["p.mci_mil"]], t = t)
  diag(t.TP2_mci) <- NA
  t.TP2_mci[1,1] <- 1-sum(t.TP2_mci[1,-1])
  t.TP2_mci[2,2] <- 1-sum(t.TP2_mci[2,-2])
  t.TP2_mci[3,3] <- 1-sum(t.TP2_mci[3,-3])
  t.TP2_mci[4,4] <- 1-sum(t.TP2_mci[4,-4])
  t.TP2_mci
  ## simple markov chain
  TP2_mci_st <- c(1,0,0,0)
  TP2_mci_st <- TP2_mci_st %*% t.TP2_mci
  TP2_mci_st <- TP2_mci_st %*% t.TP2_mci
  TP2_mci_st <- TP2_mci_st %*% t.TP2_mci
  TP2_mci_st <- TP2_mci_st %*% t.TP2_mci
  round(TP2_mci_st,3)
  
  ## handled simultaneously
  t.TP_mci <- rbind(0,cbind(0,t.TP))
  rownames(t.TP_mci)[1] <- "mci"
  colnames(t.TP_mci)[1] <- "mci"
  t.TP_mci["mci","mil"] <- l.inputs_cycle[["p.mci_mil"]]
  diag(t.TP_mci) <- NA
  t.TP_mci[1,1] <- 1-sum(t.TP_mci[1,-1])
  t.TP_mci[2,2] <- 1-sum(t.TP_mci[2,-2])
  t.TP_mci[3,3] <- 1-sum(t.TP_mci[3,-3])
  t.TP_mci[4,4] <- 1-sum(t.TP_mci[4,-4])
  t.TP_mci
  D_mci <- diag(x=eigen(t.TP_mci)$values)
  V_mci <- eigen(t.TP_mci)$vectors
  t.TP_mci <- V_mci %*% D_mci^(t) %*% solve(V_mci); t.TP_mci
  round(t.TP_mci,3)
  t.TP_mci[t.TP_mci < 0] <- 0
  t.TP_mci <- round(t.TP_mci, 3)
  colnames(t.TP_mci) <- c("mci",colnames(t.TP))
  rownames(t.TP_mci) <- c("mci",rownames(t.TP))
  rowSums(t.TP_mci)
  diag(t.TP_mci) <- NA
  t.TP_mci[1,1] <- 1-sum(t.TP_mci[1,-1])
  t.TP_mci[2,2] <- 1-sum(t.TP_mci[2,-2])
  t.TP_mci[3,3] <- 1-sum(t.TP_mci[3,-3])
  t.TP_mci[4,4] <- 1-sum(t.TP_mci[4,-4])
  t.TP_mci
  ## simple markov chain
  TP_mci_st <- c(1,0,0,0)
  TP_mci_st <- TP_mci_st %*% t.TP_mci
  TP_mci_st <- TP_mci_st %*% t.TP_mci
  TP_mci_st <- TP_mci_st %*% t.TP_mci
  TP_mci_st <- TP_mci_st %*% t.TP_mci
  round(TP_mci_st,3)
  ## compare 2 methods
  round(t.TP2_mci, 3) # TP
  round(t.TP_mci, 3)
  round(TP2_mci_st,3) # markov chain at cycle n
  round(TP_mci_st,3)
  round(TP2_mci_st - TP_mci_st,2)
  
  # adjust mortality table
  m.mortality_rate_US_2019_cycle <- m.mortality_rate_US_2019 * t # mortality table is rate, which can be multiplied (as compared to probability)
  t.rows <- rep(x = 1:nrow(m.mortality_rate_US_2019_cycle), each = 1/t) # create vector for age row numbers and repeat each age row multiple times
  m.mortality_rate_US_2019_cycle <- m.mortality_rate_US_2019_cycle[t.rows,] # select each row multiple times
  
  # adjust inputs to cycle time (list is duplicated to force check all parameters for adjustment)
  l.inputs_cycle2 <- l.inputs_cycle
  l.inputs_cycle2 <- list(
    v.names_state = c("mcion_c","mciof_c","milon_c","milof_c","mod_c","sev_c","mci_i","mil_i","mod_i","sev_i","dth"), # disease states: mci = mild cognitive impairment; mil = mild dementia; mod = moderate dementia; sev = severe dementia; dth = dead; x_i = living in institutional setting (without '_i' = living in community)
    v.names_strat = c("soc","int"), 
    age_start = l.inputs_cycle2[["age_start"]] / t, 
    sex = "weighted", 
    p.starting_state_mci = 0.55, 
    n.cycle = l.inputs_cycle2[["n.cycle"]] / t, 
    p.mci_mil = f.p_time(p = l.inputs_cycle2[["p.mci_mil"]], t = t), 
    p.mci_mod = f.p_time(p = l.inputs_cycle2[["p.mci_mod"]], t = t), 
    p.mci_sev = f.p_time(p = l.inputs_cycle2[["p.mci_sev"]], t = t), 
    p.mil_mci = 0, # assumption: changed to 0
    p.mil_mod = t.TP2["mil","mod"], 
    p.mil_sev = t.TP2["mil","sev"], 
    p.mod_mil = t.TP2["mod","mil"], 
    p.mod_sev = t.TP2["mod","sev"], 
    p.sev_mil = t.TP2["sev","mil"], 
    p.sev_mod = t.TP2["sev","mod"], 
    p.mci_i = f.p_time(p = l.inputs_cycle2[["p.mci_i"]], t = t), 
    p.mil_i = f.p_time(p = l.inputs_cycle2[["p.mil_i"]], t = t), 
    p.mod_i = f.p_time(p = l.inputs_cycle2[["p.mod_i"]], t = t), 
    p.sev_i = f.p_time(p = l.inputs_cycle2[["p.sev_i"]], t = t), 
    m.r.mortality = m.mortality_rate_US_2019_cycle, 
    hr.mort_mci = 1.82, 
    hr.mort_mil = 2.92, 
    hr.mort_mod = 3.85, 
    hr.mort_sev = 9.52, 
    rr.tx_mci_mil = 0.69, 
    rr.tx_mci_mod = 1, 
    rr.tx_mci_sev = 1, 
    rr.tx_mil_mod = 0.69, 
    rr.tx_mil_sev = 0.69, 
    rr.tx_mci_mil_dis = 1, 
    rr.tx_mci_mod_dis = 1, 
    rr.tx_mci_sev_dis = 1, 
    rr.tx_mil_mod_dis = 1, 
    rr.tx_mil_sev_dis = 1, 
    p.tx_discontinuation1 = f.p_time(p = l.inputs_cycle2[["p.tx_discontinuation1"]], t = t), 
    p.tx_discontinuation2 = f.p_time(p = l.inputs_cycle2[["p.tx_discontinuation2"]], t = t), 
    tx_discontinuation2_begin = l.inputs_cycle2[["tx_discontinuation2_begin"]] / t, 
    tx_duration = l.inputs_cycle2[["tx_duration"]] / t, 
    tx_waning = f.p_time(p = l.inputs_cycle2[["tx_waning"]], t = t), 
    tx_waning_dis = f.p_time(p = l.inputs_cycle2[["tx_waning_dis"]], t = t), 
    u.mci_pt = l.inputs_cycle2[["u.mci_pt"]] * t, 
    u.mil_pt = l.inputs_cycle2[["u.mil_pt"]] * t, 
    u.mod_pt = l.inputs_cycle2[["u.mod_pt"]] * t, 
    u.sev_pt = l.inputs_cycle2[["u.sev_pt"]] * t, 
    u.mci_pt_i = l.inputs_cycle2[["u.mci_pt_i"]] * t, 
    u.mil_pt_i = l.inputs_cycle2[["u.mil_pt_i"]] * t, 
    u.mod_pt_i = l.inputs_cycle2[["u.mod_pt_i"]] * t, 
    u.sev_pt_i = l.inputs_cycle2[["u.sev_pt_i"]] * t, 
    u.mci_ic = l.inputs_cycle2[["u.mci_ic"]] * t, 
    u.mil_ic = l.inputs_cycle2[["u.mil_ic"]] * t, 
    u.mod_ic = l.inputs_cycle2[["u.mod_ic"]] * t, 
    u.sev_ic = l.inputs_cycle2[["u.sev_ic"]] * t, 
    u.mci_ic_i = l.inputs_cycle2[["u.mci_ic_i"]] * t, 
    u.mil_ic_i = l.inputs_cycle2[["u.mil_ic_i"]] * t, 
    u.mod_ic_i = l.inputs_cycle2[["u.mod_ic_i"]] * t, 
    u.sev_ic_i = l.inputs_cycle2[["u.sev_ic_i"]] * t, 
    u.Tx_start = l.inputs_cycle2[["u.Tx_start"]], # not affected by cycle time
    c.mci_hc = l.inputs_cycle2[["c.mci_hc"]] * t, 
    c.mil_hc = l.inputs_cycle2[["c.mil_hc"]] * t, 
    c.mod_hc = l.inputs_cycle2[["c.mod_hc"]] * t, 
    c.sev_hc = l.inputs_cycle2[["c.sev_hc"]] * t, 
    c.mci_hc_i = l.inputs_cycle2[["c.mci_hc_i"]] * t, 
    c.mil_hc_i = l.inputs_cycle2[["c.mil_hc_i"]] * t, 
    c.mod_hc_i = l.inputs_cycle2[["c.mod_hc_i"]] * t, 
    c.sev_hc_i = l.inputs_cycle2[["c.sev_hc_i"]] * t, 
    c.mci_sc = l.inputs_cycle2[["c.mci_sc"]] * t, 
    c.mil_sc = l.inputs_cycle2[["c.mil_sc"]] * t, 
    c.mod_sc = l.inputs_cycle2[["c.mod_sc"]] * t, 
    c.sev_sc = l.inputs_cycle2[["c.sev_sc"]] * t, 
    c.mci_sc_i = l.inputs_cycle2[["c.mci_sc_i"]] * t, 
    c.mil_sc_i = l.inputs_cycle2[["c.mil_sc_i"]] * t, 
    c.mod_sc_i = l.inputs_cycle2[["c.mod_sc_i"]] * t, 
    c.sev_sc_i = l.inputs_cycle2[["c.sev_sc_i"]] * t, 
    c.mci_ic =  l.inputs_cycle2[["c.mci_ic"]] * t, 
    c.mil_ic = l.inputs_cycle2[["c.mil_ic"]] * t, 
    c.mod_ic = l.inputs_cycle2[["c.mod_ic"]] * t, 
    c.sev_ic = l.inputs_cycle2[["c.sev_ic"]] * t, 
    c.mci_ic_i = l.inputs_cycle2[["c.mci_ic_i"]] * t, 
    c.mil_ic_i = l.inputs_cycle2[["c.mil_ic_i"]] * t, 
    c.mod_ic_i = l.inputs_cycle2[["c.mod_ic_i"]] * t, 
    c.sev_ic_i = l.inputs_cycle2[["c.sev_ic_i"]] * t, 
    c.Tx = l.inputs_cycle2[["c.Tx"]] * t, 
    c.Tx_start = l.inputs_cycle2[["c.Tx_start"]], # not affected by cycle time
    discount_EFFECT = l.inputs_cycle2[["discount_EFFECT"]] * t, 
    discount_QALY = l.inputs_cycle2[["discount_QALY"]] * t, 
    discount_COST = l.inputs_cycle2[["discount_COST"]] * t, 
    wtp = l.inputs_cycle2[["wtp"]], # not affected by cycle time
    half_cycle_correction = FALSE # could be no longer needed with small cycle time
  )
  
  # run scenario
  l.out_cycle2 <- f.run_scenario(l.inputs = l.inputs_cycle2, detailed = TRUE)
  
  # additional results
  m.result_cycle2 <- matrix(data = NA, nrow = ncol(l.out_cycle2$l.out$soc$m.out), ncol = 3, dimnames = list(colnames(l.out_cycle2$l.out$soc$m.out), c("soc","int","dif")))
  m.result_cycle2[,"soc"] <- colSums(l.out_cycle2$l.out$soc$m.out)
  m.result_cycle2[,"int"] <- colSums(l.out_cycle2$l.out$int$m.out)
  m.result_cycle2[,"dif"] <- m.result_cycle2[,"int"] - m.result_cycle2[,"soc"]
  
  # store result in sensitivity analysis table
  m.table1[12,"mcimil"] <- sum(m.result_cycle2[c("mci","mil"),"dif"])/(1/t)
  m.table1[12,"ly"] <- m.result_cycle2["ly","dif"]/(1/t)
  m.table1[12,"qaly"] <- m.result_cycle2["qaly","dif"]
  m.table1[12,"cost_dxtx"] <- sum(m.result_cycle2[c("cost_dx","cost_tx"),"dif"])
  m.table1[12,"cost_care"] <- sum(m.result_cycle2[c("cost_hc","cost_sc","cost_ic"),"dif"])
  m.table1[12,"nhb"] <- m.result_cycle2["nhb","dif"]
  m.table1[12,"icer"] <- calculate_icers(cost = l.out_cycle2[["df.out"]][,"COST"], effect = l.out_cycle2[["df.out"]][,"QALY"], strategies = l.out_cycle2[["df.out"]][,"strategy"])[2,"ICER"]
  
}



######################################## 5.3.2. CALIBRATE TO TIME SHIFT ########################################

if(F) {
  
  # copy input estimates for calibration (cycle time must be 1/24 year, i.e., 0.5 month)
  l.inputs_cal <- l.inputs_cycle2
  
  # (temporary value for debugging; can be ignored)
  x<-0.69
  
  # function to estimate proportion at the required time delay
  f.headroom <- function(x) {
    l.inputs_cal[["rr.tx_mci_mil"]] <- x
    l.inputs_cal[["rr.tx_mil_mod"]] <- x
    l.inputs_cal[["rr.tx_mil_sev"]] <- x
    
    # run scenario
    l.out_cal <- f.run_scenario(l.inputs = l.inputs_cal, detailed=TRUE)
    
    # additional results
    m.result_cal <- matrix(data = NA, nrow = ncol(l.out_cal$l.out$soc$m.out), ncol = 3, dimnames = list(colnames(l.out_cal$l.out$soc$m.out), c("soc","int","dif")))
    m.result_cal[,"soc"] <- colSums(l.out_cal$l.out$soc$m.out)
    m.result_cal[,"int"] <- colSums(l.out_cal$l.out$int$m.out)
    m.result_cal[,"dif"] <- m.result_cal[,"int"] - m.result_cal[,"soc"]
    
    # temporary result
    round(l.out_cal[["l.out"]][["soc"]][["m.out"]][1:37,c("mci","mil","mod","sev","dth")],3) # soc: state trace up to cycle 37 (=1.5 years with 1/24 cycle length)
    round(l.out_cal[["l.out"]][["int"]][["m.out"]][1:37,c("mci","mil","mod","sev","dth")],3) # int: idem
    round(sum(l.out_cal[["l.out"]][["soc"]][["m.out"]][1:36+1,c("mci","mil")]),3)/2 # soc: person-months up to 18 months (at 1/24 cycle length)
    round(sum(l.out_cal[["l.out"]][["int"]][["m.out"]][1:36+1,c("mci","mil")]),3)/2 # int: person-months up to 18 months (at 1/24 cycle length)
    
    # pick: 
    # 1A: obtain to-be calibrated result: proportion in mci at 18 months
    if(T) {
      p.int18m <- l.out_cal[["l.out"]][["int"]][["m.out"]][36+1,"mci"] # proportion mci int at cycle 37 (=1.5 year, which matches trial duration); change mci to mil, see below
      p.soc18m <- l.out_cal[["l.out"]][["soc"]][["m.out"]][36+1-2*5.5,"mci"] # proportion mci soc at cycle 26 (=1.5 year minus 5.5 months, which matches the to-be calibrated time delay); change mci to mil, see below
      round(p.int18m,3)
      round(p.soc18m,3)
      dif <- p.int18m - p.soc18m # calculate difference
    }
    # 1B: obtain to-be calibrated result: proportion in mil at 18 months
    if(F) {
      p.int18m <- l.out_cal[["l.out"]][["int"]][["m.out"]][36+1,"mil"] # proportion mil int at cycle 37 (=1.5 year, which matches trial duration); change mci to mil, see below
      p.soc18m <- l.out_cal[["l.out"]][["soc"]][["m.out"]][36+1-2*5.5,"mil"] # proportion mil soc at cycle 26 (=1.5 year minus 5.5 months, which matches the to-be calibrated time delay); change mci to mil, see below
      round(p.int18m,3)
      round(p.soc18m,3)
      dif <- p.int18m - p.soc18m # calculate difference
    }
    
    # return difference
    return(dif)
  }
  
  # test run the function
  f.headroom(x = 0.1)
  
  # calibrate relative risk treatment effect such that proportion mci in the int strategy at 18 months is the same as proportion mci in the soc strategy at 18-5.5=12.5 months (i.e., time shift of 5.5 months)
  result_optimize <- optimize(f = function(x) abs(f.headroom(x)), interval = c(0.01,0.99), tol = 0.0001, maximum = FALSE)
  print(result_optimize) # 'minimum' gives the rr at which the time shift is 5.5 months
  # result: for mci the rr = 0.581 and for mil the rr = 0.547 (after manually adjusting mci to mil in the calibration function and rerunning the code); average = 
  mean(c(0.581,0.547)) # 0.56
  
  # run scenario at treatment rr = 0.56 (applied to icer base case with cycle length 1 year)
  l.inputs_icer3 <- l.inputs_icer
  l.inputs_icer3[["rr.tx_mci_mil"]] <- 0.56
  l.inputs_icer3[["rr.tx_mil_mod"]] <- 0.56
  l.inputs_icer3[["rr.tx_mil_sev"]] <- 0.56
  l.out_icer3 <- f.run_scenario(l.inputs = l.inputs_icer3, detailed = TRUE)
  
  # additional results
  m.result_icer3 <- matrix(data = NA, nrow = ncol(l.out_icer3$l.out$soc$m.out), ncol = 3, dimnames = list(colnames(l.out_icer3$l.out$soc$m.out), c("soc","int","dif")))
  m.result_icer3[,"soc"] <- colSums(l.out_icer3$l.out$soc$m.out)
  m.result_icer3[,"int"] <- colSums(l.out_icer3$l.out$int$m.out)
  m.result_icer3[,"dif"] <- m.result_icer3[,"int"] - m.result_icer3[,"soc"]
  
  # store result in sensitivity analysis table
  m.table1[3,"mcimil"] <- sum(m.result_icer3[c("mci","mil"),"dif"])
  m.table1[3,"ly"] <- m.result_icer3["ly","dif"]
  m.table1[3,"qaly"] <- m.result_icer3["qaly","dif"]
  m.table1[3,"cost_dxtx"] <- sum(m.result_icer3[c("cost_dx","cost_tx"),"dif"])
  m.table1[3,"cost_care"] <- sum(m.result_icer3[c("cost_hc","cost_sc","cost_ic"),"dif"])
  m.table1[3,"nhb"] <- m.result_icer3["nhb","dif"]
  m.table1[3,"icer"] <- calculate_icers(cost = l.out_icer3[["df.out"]][,"COST"], effect = l.out_icer3[["df.out"]][,"QALY"], strategies = l.out_icer3[["df.out"]][,"strategy"])[2,"ICER"]
  
  # copy sensitivity analysis results table to clipboard
    # write.table(x = round(m.table1,2), file = "clipboard", sep = "\t", row.names = FALSE, col.names = FALSE)
  
}



######################################## 5.3.3. REPLICATION: HERRING ########################################

if(F) {
  
  # U.S. general population life table 2017
  m.lifetable_US_2017 <- as.matrix(read.csv(file="life_tables/lifetable_US_2017.csv", header=TRUE))[,c("male","female","total")]
  ## convert probability to rate
  m.mortality_rate_US_2017 <- -log(1-(m.lifetable_US_2017))
  ## weight rate for male and female
  m.mortality_rate_US_2017 <- cbind(m.mortality_rate_US_2017, weighted=NA)
  m.mortality_rate_US_2017[,"weighted"] <- m.mortality_rate_US_2017[,"male"]*(1-0.524) + m.mortality_rate_US_2017[,"female"]*0.524
  
  # input parameters
  l.inputs_herring <- list(
    v.names_state = c("mcion_c","mciof_c","milon_c","milof_c","mod_c","sev_c","mci_i","mil_i","mod_i","sev_i","dth"), # disease states: mci = mild cognitive impairment; mil = mild dementia; mod = moderate dementia; sev = severe dementia; dth = dead; x_i = living in institutional setting (without '_i' = living in community)
    v.names_strat = c("soc","int"), 
    age_start = 65, 
    sex = "weighted", 
    p.starting_state_mci = 1, 
    n.cycle = 35, 
    p.mci_mil = 0.232 * 0.727, 
    p.mci_mod = 0.232 * 0.273, 
    p.mci_sev = 0, 
    p.mil_mci = 0.033, 
    p.mil_mod = 0.352, 
    p.mil_sev = 0.044, 
    p.mod_mil = 0.029, 
    p.mod_sev = 0.42, 
    p.sev_mil = 0, 
    p.sev_mod = 0.019, 
    p.mci_i = 0, 
    p.mil_i = 0, 
    p.mod_i = 0, 
    p.sev_i = 0, 
    m.r.mortality = m.mortality_rate_US_2017, 
    hr.mort_mci = 1.48, 
    hr.mort_mil = 2.84, 
    hr.mort_mod = 2.84, 
    hr.mort_sev = 2.84, 
    rr.tx_mci_mil = 0.69, 
    rr.tx_mci_mod = 0.69, 
    rr.tx_mci_sev = 0.69, 
    rr.tx_mil_mod = 0.69, 
    rr.tx_mil_sev = 0.69, 
    rr.tx_mci_mil_dis = 1, 
    rr.tx_mci_mod_dis = 1, 
    rr.tx_mci_sev_dis = 1, 
    rr.tx_mil_mod_dis = 1, 
    rr.tx_mil_sev_dis = 1, 
    p.tx_discontinuation1 = 0, 
    p.tx_discontinuation2 = 0, 
    tx_discontinuation2_begin = 29, 
    tx_duration = 29, 
    tx_waning = 0, 
    tx_waning_dis = 0, 
    u.mci_pt = 0.800, 
    u.mil_pt = 0.740, 
    u.mod_pt = 0.590, 
    u.sev_pt = 0.360, 
    u.mci_pt_i = 0, 
    u.mil_pt_i = 0, 
    u.mod_pt_i = 0, 
    u.sev_pt_i = 0, 
    u.mci_ic = 0, 
    u.mil_ic = -0.036, 
    u.mod_ic = -0.070, 
    u.sev_ic = -0.086, 
    u.mci_ic_i = 0, 
    u.mil_ic_i = 0, 
    u.mod_ic_i = 0, 
    u.sev_ic_i = 0, 
    u.Tx_start = 0, 
    c.mci_hc = 0, 
    c.mil_hc = 0, 
    c.mod_hc = 0, 
    c.sev_hc = 0, 
    c.mci_hc_i = 0, 
    c.mil_hc_i = 0, 
    c.mod_hc_i = 0, 
    c.sev_hc_i = 0, 
    c.mci_sc = 0, 
    c.mil_sc = 0, 
    c.mod_sc = 0, 
    c.sev_sc = 0, 
    c.mci_sc_i = 0, 
    c.mil_sc_i = 0, 
    c.mod_sc_i = 0, 
    c.sev_sc_i = 0, 
    c.mci_ic =  0, 
    c.mil_ic = 0, 
    c.mod_ic = 0, 
    c.sev_ic = 0, 
    c.mci_ic_i = 0, 
    c.mil_ic_i = 0, 
    c.mod_ic_i = 0, 
    c.sev_ic_i = 0, 
    c.Tx = 0, 
    c.Tx_start = 0, 
    discount_EFFECT = 0, 
    discount_QALY = 0.03, 
    discount_COST = 0, 
    wtp = 0, 
    half_cycle_correction = TRUE
  )
  
  # run scenario
  l.out_herring <- f.run_scenario(l.inputs = l.inputs_herring, detailed = TRUE)
  
  # additional results
  m.result_herring <- matrix(data = NA, nrow = ncol(l.out_herring$l.out$soc$m.out), ncol = 3, dimnames = list(colnames(l.out_herring$l.out$soc$m.out), c("soc","int","dif")))
  m.result_herring[,"soc"] <- colSums(l.out_herring$l.out$soc$m.out)
  m.result_herring[,"int"] <- colSums(l.out_herring$l.out$int$m.out)
  m.result_herring[,"dif"] <- m.result_herring[,"int"] - m.result_herring[,"soc"]
  print(round(m.result_herring[c("ly","qaly","cost"),c("soc","int","dif")],2))
  
  # compare to publication
  print(round(m.result_herring[c("ly","ontx","mci"),c("soc","int","dif")],2))
  # dem
  print(round(sum(m.result_herring[c("mil","mod","sev"),"soc"]),2))
  print(round(sum(m.result_herring[c("mil","mod","sev"),"int"]),2))
  print(round(sum(m.result_herring[c("mil","mod","sev"),"int"]) - sum(m.result_herring[c("mil","mod","sev"),"soc"]),2))
  
}


######################################## 5.4. CONFERENCE ALZHEIMER EUROPE 2024 ########################################

if(F) {
  
  # copy inputs ICER replication
  l.inputs_eu <- l.inputs_icer
  
  # half-cycle correction
  l.inputs_eu[["half_cycle_correction"]] <- FALSE
  
  # disease progression MCI to dementia in Amyloid positive & neuronal loss undetermined [Vos, 2015: https://doi.org/10.1093/brain/awv029] 
  ## Operationalized by diagnostic criteria NIA-AA categories: 'NIA-AA high AD' (Amyloid+, Injury+) and 'conflicting IAP' (Amyloid+, Injury-)
  ## corresponding 3-year cumulative incidence probability: 'high AD' = 59% (AD dementia; table 3) and 4% (non-AD dementia; mentioned in text), and 22% (AD dementia; table 3) and 4% (non-AD dementia; mentioned in text) with prevalence of 353 and 49 respectively (respectively)
  ## This results into a weighted 3-year cumulative incidence of (converting all 4 probabilities to rates before weighting and averaging):
  temp.est2 <- 1-exp(- ( (-log(1-0.59) + -log(1-0.04))*353 + (-log(1-0.22) + -log(1-0.04))*49 ) / (353+49) )
  ## and corresponding 1-year probability of 
  temp.est2 <- 1-exp(- -log(1-temp.est2)/3)
  temp.est2
  l.inputs_eu[["p.mci_mil"]] <- round(temp.est2,3) # 0.248
  l.inputs_eu[["p.mil_mci"]] <- 0
  
  # disease progression dementia [Wimo, 2020: https://doi.org/10.3233/jad-191055]
  l.inputs_eu[["p.mci_mod"]] <- 0
  l.inputs_eu[["p.mci_sev"]] <- 0
  l.inputs_eu[["p.mil_mci"]] <- 0
  l.inputs_eu[["p.mil_mod"]] <- 0.293
  l.inputs_eu[["p.mil_sev"]] <- 0.001
  l.inputs_eu[["p.mod_mil"]] <- 0.087
  l.inputs_eu[["p.mod_sev"]] <- 0.109
  l.inputs_eu[["p.sev_mil"]] <- 0.000
  l.inputs_eu[["p.sev_mod"]] <- 0.196
  
  # institutionalization set to 0 (because rate not obtained for EU regions/countries)
  l.inputs_eu[["p.mci_i"]] <- 0
  l.inputs_eu[["p.mil_i"]] <- 0
  l.inputs_eu[["p.mod_i"]] <- 0
  l.inputs_eu[["p.sev_i"]] <- 0
  
  # mortality [Wimo, 2020: https://doi.org/10.3233/jad-191055]
  l.inputs_eu[["hr.mort_mci"]] <- 1
  l.inputs_eu[["hr.mort_mil"]] <- 1.318 * 1.82
  l.inputs_eu[["hr.mort_mod"]] <- 2.419 * 1.82
  l.inputs_eu[["hr.mort_sev"]] <- 4.267 * 1.82
  
  # utilities
  l.inputs_eu[["u.mci_pt"]] <- mean(c(0.7 ,0.75)) # Landeiro review 2020 studies Jonsson 2006 and Hessman 2016 (as these are EU mixed setting including MCI and dementia stages)
  l.inputs_eu[["u.mil_pt"]] <- mean(c(0.65,0.61))
  l.inputs_eu[["u.mod_pt"]] <- mean(c(0.51,0.41))
  l.inputs_eu[["u.sev_pt"]] <- mean(c(0.40,0.21))
  l.inputs_eu[["u.mci_pt_i"]] <- 0
  l.inputs_eu[["u.mil_pt_i"]] <- 0
  l.inputs_eu[["u.mod_pt_i"]] <- 0
  l.inputs_eu[["u.sev_pt_i"]] <- 0
  l.inputs_eu[["u.mci_ic"]] <- 0
  l.inputs_eu[["u.mil_ic"]] <- 0
  l.inputs_eu[["u.mod_ic"]] <- 0
  l.inputs_eu[["u.sev_ic"]] <- 0
  l.inputs_eu[["u.mci_ic_i"]] <- 0
  l.inputs_eu[["u.mil_ic_i"]] <- 0
  l.inputs_eu[["u.mod_ic_i"]] <- 0
  l.inputs_eu[["u.sev_ic_i"]] <- 0
  
  # set costs to 0
  l.inputs_eu[["c.mci_hc"]] <- 0
  l.inputs_eu[["c.mil_hc"]] <- 0
  l.inputs_eu[["c.mod_hc"]] <- 0
  l.inputs_eu[["c.sev_hc"]] <- 0
  l.inputs_eu[["c.mci_hc_i"]] <- 0
  l.inputs_eu[["c.mil_hc_i"]] <- 0
  l.inputs_eu[["c.mod_hc_i"]] <- 0
  l.inputs_eu[["c.sev_hc_i"]] <- 0
  l.inputs_eu[["c.mci_sc"]] <- 0 
  l.inputs_eu[["c.mil_sc"]] <- 0 
  l.inputs_eu[["c.mod_sc"]] <- 0 
  l.inputs_eu[["c.sev_sc"]] <- 0 
  l.inputs_eu[["c.mci_sc_i"]] <- 0
  l.inputs_eu[["c.mil_sc_i"]] <- 0
  l.inputs_eu[["c.mod_sc_i"]] <- 0
  l.inputs_eu[["c.sev_sc_i"]] <- 0
  l.inputs_eu[["c.mci_ic"]] <-  0
  l.inputs_eu[["c.mil_ic"]] <- 0
  l.inputs_eu[["c.mod_ic"]] <- 0
  l.inputs_eu[["c.sev_ic"]] <- 0
  l.inputs_eu[["c.mci_ic_i"]] <- 0
  l.inputs_eu[["c.mil_ic_i"]] <- 0
  l.inputs_eu[["c.mod_ic_i"]] <- 0
  l.inputs_eu[["c.sev_ic_i"]] <- 0
  l.inputs_eu[["c.Tx"]] <- (26500 + (52/2)*78.35) / 1.1827 # treatment costs (1.1827 exchange rate Euro to US dollar from https://ec.europa.eu/eurostat/databrowser/view/tec00033/default/table?lang=en&category=t_ert)
  l.inputs_eu[["c.Tx_start"]] <- (261.10*4 + 261.10*3*0.215) / 1.1827
  
  # load life tables from selection of EU countries
  a.lifetable <- array(data=NA, dim=c(100,3,5), dimnames=list(NULL,c("male","female","weighted"),c("ES","NL","PL","SE","UK")))
  a.lifetable[,c("male","female"),"ES"] <- as.matrix(read.csv(file="life_tables/lifetable_ES_2021.csv", header=TRUE))[,c("male","female")]
  a.lifetable[,c("male","female"),"NL"] <- as.matrix(read.csv(file="life_tables/lifetable_NL_2021.csv", header=TRUE))[,c("male","female")]
  a.lifetable[,c("male","female"),"PL"] <- as.matrix(read.csv(file="life_tables/lifetable_PL_2021.csv", header=TRUE))[,c("male","female")]
  a.lifetable[,c("male","female"),"SE"] <- as.matrix(read.csv(file="life_tables/lifetable_SE_2021.csv", header=TRUE))[,c("male","female")]
  a.lifetable[,c("male","female"),"UK"] <- as.matrix(read.csv(file="life_tables/lifetable_UK_2021.csv", header=TRUE))[,c("male","female")]
  matplot(x=a.lifetable[1:99,"male",], type="l", col=rainbow(5))
  legend(x="topleft", legend=c("ES","NL","PL","SE","UK"), col=rainbow(5), lty=c(1:5))
  a.lifetable[,"weighted",] <- a.lifetable[,"male",] * 0.48 + a.lifetable[,"female",] * 0.52 # weights (same as ICER replication)
  
  # region-specific inputs: life table
  l.inputs_eas <- l.inputs_nor <- l.inputs_wes <- l.inputs_sou <- l.inputs_bri <- l.inputs_eu
  l.inputs_eas[["m.r.mortality"]] <- a.lifetable[,,"PL"]
  l.inputs_nor[["m.r.mortality"]] <- a.lifetable[,,"SE"]
  l.inputs_wes[["m.r.mortality"]] <- a.lifetable[,,"NL"]
  l.inputs_sou[["m.r.mortality"]] <- a.lifetable[,,"ES"]
  l.inputs_bri[["m.r.mortality"]] <- a.lifetable[,,"UK"]
  
  # costs in EU regions [Jonsson, 2023: https://doi.org/10.1007/s40273-022-01212-z supplemental material]
  c_na <- c(NA, NA, NA, NA, NA, NA, NA)
  c_e_mil <- c(125.69, 389.31, 819.81, 1952.78, 1023.69, 78.33, 3226.81)
  c_e_mod <- c(371.75, 437.27, 851.49, 2627.22, 1870.16, 478.9, 3032.98)
  c_e_sev <- c(748.95, 533.84, 1068.32, 4611.37, 857.21, 665.6, 2750.78)
  c_n_mil <- c(2247.4, 986.56, 2007.24, 6075.43, 4440.07, 504.38, 4614.44)
  c_n_mod <- c(2121.17, 1061.47, 1288.69, 21130.93, 6107.78, 717.71, 5112.33)
  c_n_sev <- c(1904.31, 1005.6, 1130.67, 40880.29, 7175.8, 305.59, 5795.96)
  c_w_mil <- c(1806.92, 1371.11, 1017.85, 4449.83, 8098.58, 988.41, 14251.28)
  c_w_mod <- c(2175.48, 2355.42, 1233.83, 13708.92, 8379.1, 1135.68, 18946)
  c_w_sev <- c(1672.04, 2997.69, 1154.24, 15763.69, 9191.27, 654.06, 24671.03)
  c_s_mil <- c(326.41, 870.02, 1584.38, 133.05, 2433.26, 465.42, 14607.64)
  c_s_mod <- c(901.69, 1279.84, 1748.57, 959.29, 4847.88, 1632.47, 29582.94)
  c_s_sev <- c(1743.6, 1158.02, 2113.16, 7617.24, 4541.32, 638.07, 44094.34)
  c_b_mil <- c(4048.92, 889.22, 728.12, 2562.5, 2446.38, 541.54, 8692.37)
  c_b_mod <- c(2225.86, 1373.19, 593.19, 6909.16, 4661.5, 725.03, 17735.22)
  c_b_sev <- c(2277.26, 1074.67, 967.02, 15490.36, 5687.87, 2088.04, 34372.29)
  a.c_eu <- array(
    data = c(
      c_na, c_e_mil, c_e_mod, c_e_sev, 
      c_na, c_n_mil, c_n_mod, c_n_sev, 
      c_na, c_w_mil, c_w_mod, c_w_sev, 
      c_na, c_s_mil, c_s_mod, c_s_sev, 
      c_na, c_b_mil, c_b_mod, c_b_sev
    ), 
    dim = c(7,4,5), 
    dimnames = list(c("inp","out","pha","ins","com","icw","ict"),c("mci","mil","mod","sev"),c("eas","nor","wes","sou","bri"))
  ); a.c_eu
  barplot(height = a.c_eu[,,"wes"]) # check with original publication figure 3
  a.c_eu[,"mci",] <- a.c_eu[,"mil",] * (1.12/1.56) # add costs for MCI (assumed ratio health care sector costs between m.mil and m.mci from ICER replication)
  round(a.c_eu[c("inp","out","pha","ins","com"),,],0)
  
  # region-specific inputs: costs
  l.inputs_eas[["c.mci_hc"]] <- sum(a.c_eu[c("inp","out","pha"),"mci","eas"])
  l.inputs_eas[["c.mil_hc"]] <- sum(a.c_eu[c("inp","out","pha"),"mil","eas"])
  l.inputs_eas[["c.mod_hc"]] <- sum(a.c_eu[c("inp","out","pha"),"mod","eas"])
  l.inputs_eas[["c.sev_hc"]] <- sum(a.c_eu[c("inp","out","pha"),"sev","eas"])
  l.inputs_eas[["c.mci_sc"]] <- sum(a.c_eu[c("ins","com"),"mci","eas"])
  l.inputs_eas[["c.mil_sc"]] <- sum(a.c_eu[c("ins","com"),"mil","eas"])
  l.inputs_eas[["c.mod_sc"]] <- sum(a.c_eu[c("ins","com"),"mod","eas"])
  l.inputs_eas[["c.sev_sc"]] <- sum(a.c_eu[c("ins","com"),"sev","eas"])
  l.inputs_eas[["c.mci_ic"]] <- sum(a.c_eu[c("ict"),"mci","eas"])
  l.inputs_eas[["c.mil_ic"]] <- sum(a.c_eu[c("ict"),"mil","eas"])
  l.inputs_eas[["c.mod_ic"]] <- sum(a.c_eu[c("ict"),"mod","eas"])
  l.inputs_eas[["c.sev_ic"]] <- sum(a.c_eu[c("ict"),"sev","eas"])
  l.inputs_nor[["c.mci_hc"]] <- sum(a.c_eu[c("inp","out","pha"),"mci","nor"])
  l.inputs_nor[["c.mil_hc"]] <- sum(a.c_eu[c("inp","out","pha"),"mil","nor"])
  l.inputs_nor[["c.mod_hc"]] <- sum(a.c_eu[c("inp","out","pha"),"mod","nor"])
  l.inputs_nor[["c.sev_hc"]] <- sum(a.c_eu[c("inp","out","pha"),"sev","nor"])
  l.inputs_nor[["c.mci_sc"]] <- sum(a.c_eu[c("ins","com"),"mci","nor"])
  l.inputs_nor[["c.mil_sc"]] <- sum(a.c_eu[c("ins","com"),"mil","nor"])
  l.inputs_nor[["c.mod_sc"]] <- sum(a.c_eu[c("ins","com"),"mod","nor"])
  l.inputs_nor[["c.sev_sc"]] <- sum(a.c_eu[c("ins","com"),"sev","nor"])
  l.inputs_nor[["c.mci_ic"]] <- sum(a.c_eu[c("ict"),"mci","nor"])
  l.inputs_nor[["c.mil_ic"]] <- sum(a.c_eu[c("ict"),"mil","nor"])
  l.inputs_nor[["c.mod_ic"]] <- sum(a.c_eu[c("ict"),"mod","nor"])
  l.inputs_nor[["c.sev_ic"]] <- sum(a.c_eu[c("ict"),"sev","nor"])
  l.inputs_wes[["c.mci_hc"]] <- sum(a.c_eu[c("inp","out","pha"),"mci","wes"])
  l.inputs_wes[["c.mil_hc"]] <- sum(a.c_eu[c("inp","out","pha"),"mil","wes"])
  l.inputs_wes[["c.mod_hc"]] <- sum(a.c_eu[c("inp","out","pha"),"mod","wes"])
  l.inputs_wes[["c.sev_hc"]] <- sum(a.c_eu[c("inp","out","pha"),"sev","wes"])
  l.inputs_wes[["c.mci_sc"]] <- sum(a.c_eu[c("ins","com"),"mci","wes"])
  l.inputs_wes[["c.mil_sc"]] <- sum(a.c_eu[c("ins","com"),"mil","wes"])
  l.inputs_wes[["c.mod_sc"]] <- sum(a.c_eu[c("ins","com"),"mod","wes"])
  l.inputs_wes[["c.sev_sc"]] <- sum(a.c_eu[c("ins","com"),"sev","wes"])
  l.inputs_wes[["c.mci_ic"]] <- sum(a.c_eu[c("ict"),"mci","wes"])
  l.inputs_wes[["c.mil_ic"]] <- sum(a.c_eu[c("ict"),"mil","wes"])
  l.inputs_wes[["c.mod_ic"]] <- sum(a.c_eu[c("ict"),"mod","wes"])
  l.inputs_wes[["c.sev_ic"]] <- sum(a.c_eu[c("ict"),"sev","wes"])
  l.inputs_sou[["c.mci_hc"]] <- sum(a.c_eu[c("inp","out","pha"),"mci","sou"])
  l.inputs_sou[["c.mil_hc"]] <- sum(a.c_eu[c("inp","out","pha"),"mil","sou"])
  l.inputs_sou[["c.mod_hc"]] <- sum(a.c_eu[c("inp","out","pha"),"mod","sou"])
  l.inputs_sou[["c.sev_hc"]] <- sum(a.c_eu[c("inp","out","pha"),"sev","sou"])
  l.inputs_sou[["c.mci_sc"]] <- sum(a.c_eu[c("ins","com"),"mci","sou"])
  l.inputs_sou[["c.mil_sc"]] <- sum(a.c_eu[c("ins","com"),"mil","sou"])
  l.inputs_sou[["c.mod_sc"]] <- sum(a.c_eu[c("ins","com"),"mod","sou"])
  l.inputs_sou[["c.sev_sc"]] <- sum(a.c_eu[c("ins","com"),"sev","sou"])
  l.inputs_sou[["c.mci_ic"]] <- sum(a.c_eu[c("ict"),"mci","sou"])
  l.inputs_sou[["c.mil_ic"]] <- sum(a.c_eu[c("ict"),"mil","sou"])
  l.inputs_sou[["c.mod_ic"]] <- sum(a.c_eu[c("ict"),"mod","sou"])
  l.inputs_sou[["c.sev_ic"]] <- sum(a.c_eu[c("ict"),"sev","sou"])
  l.inputs_bri[["c.mci_hc"]] <- sum(a.c_eu[c("inp","out","pha"),"mci","bri"])
  l.inputs_bri[["c.mil_hc"]] <- sum(a.c_eu[c("inp","out","pha"),"mil","bri"])
  l.inputs_bri[["c.mod_hc"]] <- sum(a.c_eu[c("inp","out","pha"),"mod","bri"])
  l.inputs_bri[["c.sev_hc"]] <- sum(a.c_eu[c("inp","out","pha"),"sev","bri"])
  l.inputs_bri[["c.mci_sc"]] <- sum(a.c_eu[c("ins","com"),"mci","bri"])
  l.inputs_bri[["c.mil_sc"]] <- sum(a.c_eu[c("ins","com"),"mil","bri"])
  l.inputs_bri[["c.mod_sc"]] <- sum(a.c_eu[c("ins","com"),"mod","bri"])
  l.inputs_bri[["c.sev_sc"]] <- sum(a.c_eu[c("ins","com"),"sev","bri"])
  l.inputs_bri[["c.mci_ic"]] <- sum(a.c_eu[c("ict"),"mci","bri"])
  l.inputs_bri[["c.mil_ic"]] <- sum(a.c_eu[c("ict"),"mil","bri"])
  l.inputs_bri[["c.mod_ic"]] <- sum(a.c_eu[c("ict"),"mod","bri"])
  l.inputs_bri[["c.sev_ic"]] <- sum(a.c_eu[c("ict"),"sev","bri"])
  
  # run scenarios
  l.out_eas <- f.run_scenario(l.inputs = l.inputs_eas, detailed = TRUE)
  l.out_nor <- f.run_scenario(l.inputs = l.inputs_nor, detailed = TRUE)
  l.out_wes <- f.run_scenario(l.inputs = l.inputs_wes, detailed = TRUE)
  l.out_sou <- f.run_scenario(l.inputs = l.inputs_sou, detailed = TRUE)
  l.out_bri <- f.run_scenario(l.inputs = l.inputs_bri, detailed = TRUE)
  
  # additional results
  m.result_eas <- matrix(data = NA, nrow = ncol(l.out_eas$l.out$soc$m.out), ncol = 3, dimnames = list(colnames(l.out_eas$l.out$soc$m.out), c("soc","int","dif")))
  m.result_eas[,"soc"] <- colSums(l.out_eas$l.out$soc$m.out)
  m.result_eas[,"int"] <- colSums(l.out_eas$l.out$int$m.out)
  m.result_eas <- rbind(m.result_eas, qaly_pt_mb=m.result_eas["qaly_pt",] * l.inputs_eas[["wtp"]])
  m.result_nor <- matrix(data = NA, nrow = ncol(l.out_nor$l.out$soc$m.out), ncol = 3, dimnames = list(colnames(l.out_nor$l.out$soc$m.out), c("soc","int","dif")))
  m.result_nor[,"soc"] <- colSums(l.out_nor$l.out$soc$m.out)
  m.result_nor[,"int"] <- colSums(l.out_nor$l.out$int$m.out)
  m.result_nor <- rbind(m.result_nor, qaly_pt_mb=m.result_nor["qaly_pt",] * l.inputs_nor[["wtp"]])
  m.result_wes <- matrix(data = NA, nrow = ncol(l.out_wes$l.out$soc$m.out), ncol = 3, dimnames = list(colnames(l.out_wes$l.out$soc$m.out), c("soc","int","dif")))
  m.result_wes[,"soc"] <- colSums(l.out_wes$l.out$soc$m.out)
  m.result_wes[,"int"] <- colSums(l.out_wes$l.out$int$m.out)
  m.result_wes <- rbind(m.result_wes, qaly_pt_mb=m.result_wes["qaly_pt",] * l.inputs_wes[["wtp"]])
  m.result_sou <- matrix(data = NA, nrow = ncol(l.out_sou$l.out$soc$m.out), ncol = 3, dimnames = list(colnames(l.out_sou$l.out$soc$m.out), c("soc","int","dif")))
  m.result_sou[,"soc"] <- colSums(l.out_sou$l.out$soc$m.out)
  m.result_sou[,"int"] <- colSums(l.out_sou$l.out$int$m.out)
  m.result_sou <- rbind(m.result_sou, qaly_pt_mb=m.result_sou["qaly_pt",] * l.inputs_sou[["wtp"]])
  m.result_bri <- matrix(data = NA, nrow = ncol(l.out_bri$l.out$soc$m.out), ncol = 3, dimnames = list(colnames(l.out_bri$l.out$soc$m.out), c("soc","int","dif")))
  m.result_bri[,"soc"] <- colSums(l.out_bri$l.out$soc$m.out)
  m.result_bri[,"int"] <- colSums(l.out_bri$l.out$int$m.out)
  m.result_bri <- rbind(m.result_bri, qaly_pt_mb=m.result_bri["qaly_pt",] * l.inputs_bri[["wtp"]])
  ## combine results
  a.result <- array(data = c(m.result_eas,m.result_nor,m.result_wes,m.result_sou,m.result_bri), dim = c(dim(m.result_eas),5), dimnames = list(rownames(m.result_eas),colnames(m.result_eas),c("eas","nor","wes","sou","bri")))
  ## additional results
  a.result[,"dif",] <- a.result[,"int",] - a.result[,"soc",]
  
  # RESULT: time in state
  a.result[c("mci","mil","mod","sev"),,]
  tis_eas <- apply(X = a.result[c("mci","mil","mod","sev"),c("soc","int"),"eas"], MARGIN = 2, FUN = cumsum)
  tis_nor <- apply(X = a.result[c("mci","mil","mod","sev"),c("soc","int"),"nor"], MARGIN = 2, FUN = cumsum)
  tis_wes <- apply(X = a.result[c("mci","mil","mod","sev"),c("soc","int"),"wes"], MARGIN = 2, FUN = cumsum)
  tis_sou <- apply(X = a.result[c("mci","mil","mod","sev"),c("soc","int"),"sou"], MARGIN = 2, FUN = cumsum)
  tis_bri <- apply(X = a.result[c("mci","mil","mod","sev"),c("soc","int"),"bri"], MARGIN = 2, FUN = cumsum)
  tis <- array(data = c(tis_eas, tis_nor, tis_wes, tis_sou, tis_bri), dim = c(4,2,5), dimnames = list(c("mci","mil","mod","sev"),c("soc","int"),c("east","north","west","south","brit")))
  barplot(height = tis[4,,], col = c("red","red")      , density = c(NA, NA), beside = TRUE, ylab = "person-years", main = "time in state")
  barplot(height = tis[3,,], col = c("orange","orange"), density = c(NA, NA), beside = TRUE, add = TRUE, xaxt = "n", yaxt = "n")
  barplot(height = tis[2,,], col = c("yellow","yellow"), density = c(NA, NA), beside = TRUE, add = TRUE, xaxt = "n", yaxt = "n")
  barplot(height = tis[1,,], col = c("green","green")  , density = c(NA, NA), beside = TRUE, add = TRUE, xaxt = "n", yaxt = "n")
  barplot(height = tis[4,,], col = c("black","black")  , density = c(0 , 30), beside = TRUE, add = TRUE, xaxt = "n", yaxt = "n")
  legend(x = "bottom", legend = c("mci","mil","mod","sev"), fill = c("green","yellow","orange","red"), ncol = 4, bg="white", inset = c(0, -0.43))
  
  # RESULT: QALY gain
  round(a.result["qaly","dif",],2)
  qaly <- a.result["qaly",c("soc","int"),]
  dimnames(qaly)[[2]] <- c("east","north","west","south","brit")
  barplot(height = qaly, col = c("cornflowerblue","cornflowerblue"), density = c(NA, NA), beside = TRUE, ylab = "quality-adjusted life years", main = "quality-adjusted life years")
  barplot(height = qaly, col = c("black","black"), density = c(0, 30), beside = TRUE, add = TRUE, xaxt = "n", yaxt = "n")
  
  # RESULT: range care sector savings
  a.result[c("cost_hc","cost_sc"),"dif",]
  round(colSums(a.result[c("cost_hc","cost_sc"),"dif",]),-2)
    # 144 to -7748
  
  # RESULT: diagnostic and drug cost
  a.result[c("cost_dx","cost_tx"),"dif",]
  round(colSums(a.result[c("cost_dx","cost_tx"),"dif",]),-3)
    # 113,000-128,000
  
  cost_eas <- apply(X = a.result[c("cost_hc","cost_sc","cost_dx","cost_tx"),c("soc","int"),"eas"], MARGIN = 2, FUN = cumsum)
  cost_nor <- apply(X = a.result[c("cost_hc","cost_sc","cost_dx","cost_tx"),c("soc","int"),"nor"], MARGIN = 2, FUN = cumsum)
  cost_wes <- apply(X = a.result[c("cost_hc","cost_sc","cost_dx","cost_tx"),c("soc","int"),"wes"], MARGIN = 2, FUN = cumsum)
  cost_sou <- apply(X = a.result[c("cost_hc","cost_sc","cost_dx","cost_tx"),c("soc","int"),"sou"], MARGIN = 2, FUN = cumsum)
  cost_bri <- apply(X = a.result[c("cost_hc","cost_sc","cost_dx","cost_tx"),c("soc","int"),"bri"], MARGIN = 2, FUN = cumsum)
  cost <- array(data = c(cost_eas, cost_nor, cost_wes, cost_sou, cost_bri), dim = c(4,2,5), dimnames = list(c("health","social","diagnostic","drug"),c("soc","int"),c("east","north","west","south","brit")))
  barplot(height = cost[4,,], col = c("yellow","yellow")      , density = c(NA, NA), beside = TRUE, ylab = "Euro", main = "costs accross care sectors")
  barplot(height = cost[3,,], col = c("red","red"), density = c(NA, NA), beside = TRUE, add = TRUE, xaxt = "n", yaxt = "n")
  barplot(height = cost[2,,], col = c("cornflowerblue","cornflowerblue"), density = c(NA, NA), beside = TRUE, add = TRUE, xaxt = "n", yaxt = "n")
  barplot(height = cost[1,,], col = c("purple","purple")  , density = c(NA, NA), beside = TRUE, add = TRUE, xaxt = "n", yaxt = "n")
  barplot(height = cost[4,,], col = c("black","black")  , density = c(0 , 30), beside = TRUE, add = TRUE, xaxt = "n", yaxt = "n")
  par(xpd=TRUE)
  legend(x = "bottom", legend = c("health","social","diagnostic","drug"), fill = c("purple","cornflowerblue","red","yellow"), ncol = 4, bg="white", inset = c(0, -0.43))
  
  # RESULT: ICER
  calculate_icers(cost = l.out_eas[["df.out"]][,"COST"], effect = l.out_eas[["df.out"]][,"QALY"], strategies = l.out_eas[["df.out"]][,"strategy"])
  calculate_icers(cost = l.out_nor[["df.out"]][,"COST"], effect = l.out_nor[["df.out"]][,"QALY"], strategies = l.out_nor[["df.out"]][,"strategy"])
  calculate_icers(cost = l.out_wes[["df.out"]][,"COST"], effect = l.out_wes[["df.out"]][,"QALY"], strategies = l.out_wes[["df.out"]][,"strategy"])
  calculate_icers(cost = l.out_sou[["df.out"]][,"COST"], effect = l.out_sou[["df.out"]][,"QALY"], strategies = l.out_sou[["df.out"]][,"strategy"])
  calculate_icers(cost = l.out_bri[["df.out"]][,"COST"], effect = l.out_bri[["df.out"]][,"QALY"], strategies = l.out_bri[["df.out"]][,"strategy"])
    # range: 309,000-346,000
  
  # RESULT: % sector savings
  a.result[c("cost_hc","cost_sc","cost_ic"),"dif",] / colSums(a.result[c("cost_hc","cost_sc","cost_ic"),"dif",])
  
  # headroom function
  f.headroom <- function(x, l.inputs, parameter) {
    l.inputs[[parameter]] <- x
    out <- f.run_scenario(l.inputs=l.inputs, detailed=FALSE)
    iNHB <- out[2,"NHB"] - out[1,"NHB"]
    return(iNHB)
  }
  
  # run headroom
  headroom_eas <- optimize(f=function(x) abs(f.headroom(x, l.inputs=l.inputs_eas, parameter="c.Tx")), interval=c(1,100000))[["minimum"]]; headroom_eas
  headroom_nor <- optimize(f=function(x) abs(f.headroom(x, l.inputs=l.inputs_nor, parameter="c.Tx")), interval=c(1,100000))[["minimum"]]; headroom_nor
  headroom_wes <- optimize(f=function(x) abs(f.headroom(x, l.inputs=l.inputs_wes, parameter="c.Tx")), interval=c(1,100000))[["minimum"]]; headroom_wes
  headroom_sou <- optimize(f=function(x) abs(f.headroom(x, l.inputs=l.inputs_sou, parameter="c.Tx")), interval=c(1,100000))[["minimum"]]; headroom_sou
  headroom_bri <- optimize(f=function(x) abs(f.headroom(x, l.inputs=l.inputs_bri, parameter="c.Tx")), interval=c(1,100000))[["minimum"]]; headroom_bri
  # temporary for testing
  l.inputs_test <- l.inputs_eas
  l.inputs_test$c.Tx <- 8000
  l.out_test <- f.run_scenario(l.inputs = l.inputs_test, detailed = TRUE)
  l.out_test$df.out
  
  # RESULT: headroom
  round(c(headroom_eas,headroom_nor,headroom_wes,headroom_sou,headroom_bri),-2)
  
  # sensitivity analysis 1: stop 18 months without sustained effect
  l.inputs_eas_sa1 <- l.inputs_eas
  l.inputs_nor_sa1 <- l.inputs_nor
  l.inputs_wes_sa1 <- l.inputs_wes
  l.inputs_sou_sa1 <- l.inputs_sou
  l.inputs_bri_sa1 <- l.inputs_bri
  l.inputs_eas_sa1[["p.tx_discontinuation2"]] <- 1
  l.inputs_eas_sa1[["tx_discontinuation2_begin"]] <- 2 # 1=no longer treated from cycle 2 onward (thus only treated for cycle 1, which is starting population for 1 year)
  l.inputs_nor_sa1[["p.tx_discontinuation2"]] <- 1
  l.inputs_nor_sa1[["tx_discontinuation2_begin"]] <- 2
  l.inputs_wes_sa1[["p.tx_discontinuation2"]] <- 1
  l.inputs_wes_sa1[["tx_discontinuation2_begin"]] <- 2
  l.inputs_sou_sa1[["p.tx_discontinuation2"]] <- 1
  l.inputs_sou_sa1[["tx_discontinuation2_begin"]] <- 2
  l.inputs_bri_sa1[["p.tx_discontinuation2"]] <- 1
  l.inputs_bri_sa1[["tx_discontinuation2_begin"]] <- 2
  # l.inputs_eas_sa1$c.Tx <- 10000 # temporary for testing
  # l.inputs_eas_sa1$discount_COST <- 0 # temporary for testing
  # l.out_eas_sa1 <- f.run_scenario(l.inputs = l.inputs_eas_sa1, detailed = TRUE)
  # print(round(l.out_eas_sa1$l.out$int$m.trace,2))
  # print(round(l.out_eas_sa1$l.out$int$m.out,2))
  l.out_eas_sa1 <- f.run_scenario(l.inputs = l.inputs_eas_sa1, detailed = TRUE)
  l.out_nor_sa1 <- f.run_scenario(l.inputs = l.inputs_nor_sa1, detailed = TRUE)
  l.out_wes_sa1 <- f.run_scenario(l.inputs = l.inputs_wes_sa1, detailed = TRUE)
  l.out_sou_sa1 <- f.run_scenario(l.inputs = l.inputs_sou_sa1, detailed = TRUE)
  l.out_bri_sa1 <- f.run_scenario(l.inputs = l.inputs_bri_sa1, detailed = TRUE)
  calculate_icers(cost = l.out_eas_sa1[["df.out"]][,"COST"], effect = l.out_eas_sa1[["df.out"]][,"QALY"], strategies = l.out_eas_sa1[["df.out"]][,"strategy"])
  calculate_icers(cost = l.out_nor_sa1[["df.out"]][,"COST"], effect = l.out_nor_sa1[["df.out"]][,"QALY"], strategies = l.out_nor_sa1[["df.out"]][,"strategy"])
  calculate_icers(cost = l.out_wes_sa1[["df.out"]][,"COST"], effect = l.out_wes_sa1[["df.out"]][,"QALY"], strategies = l.out_wes_sa1[["df.out"]][,"strategy"])
  calculate_icers(cost = l.out_sou_sa1[["df.out"]][,"COST"], effect = l.out_sou_sa1[["df.out"]][,"QALY"], strategies = l.out_sou_sa1[["df.out"]][,"strategy"])
  calculate_icers(cost = l.out_bri_sa1[["df.out"]][,"COST"], effect = l.out_bri_sa1[["df.out"]][,"QALY"], strategies = l.out_bri_sa1[["df.out"]][,"strategy"])
    # range: 290,000-325,000
  
  # sensitivity analysis 2: stop at 18 months with sustained effect
  l.inputs_eas_sa2 <- l.inputs_eas_sa1
  l.inputs_nor_sa2 <- l.inputs_nor_sa1
  l.inputs_wes_sa2 <- l.inputs_wes_sa1
  l.inputs_sou_sa2 <- l.inputs_sou_sa1
  l.inputs_bri_sa2 <- l.inputs_bri_sa1
  l.inputs_eas_sa2[["rr.tx_mci_mil_dis"]] <- 0.69
  l.inputs_eas_sa2[["rr.tx_mil_mod_dis"]] <- 0.69
  l.inputs_eas_sa2[["rr.tx_mil_sev_dis"]] <- 0.69
  l.inputs_nor_sa2[["rr.tx_mci_mil_dis"]] <- 0.69
  l.inputs_nor_sa2[["rr.tx_mil_mod_dis"]] <- 0.69
  l.inputs_nor_sa2[["rr.tx_mil_sev_dis"]] <- 0.69
  l.inputs_wes_sa2[["rr.tx_mci_mil_dis"]] <- 0.69
  l.inputs_wes_sa2[["rr.tx_mil_mod_dis"]] <- 0.69
  l.inputs_wes_sa2[["rr.tx_mil_sev_dis"]] <- 0.69
  l.inputs_sou_sa2[["rr.tx_mci_mil_dis"]] <- 0.69
  l.inputs_sou_sa2[["rr.tx_mil_mod_dis"]] <- 0.69
  l.inputs_sou_sa2[["rr.tx_mil_sev_dis"]] <- 0.69
  l.inputs_bri_sa2[["rr.tx_mci_mil_dis"]] <- 0.69
  l.inputs_bri_sa2[["rr.tx_mil_mod_dis"]] <- 0.69
  l.inputs_bri_sa2[["rr.tx_mil_sev_dis"]] <- 0.69
  l.out_eas_sa2 <- f.run_scenario(l.inputs = l.inputs_eas_sa2, detailed = TRUE)
  l.out_nor_sa2 <- f.run_scenario(l.inputs = l.inputs_nor_sa2, detailed = TRUE)
  l.out_wes_sa2 <- f.run_scenario(l.inputs = l.inputs_wes_sa2, detailed = TRUE)
  l.out_sou_sa2 <- f.run_scenario(l.inputs = l.inputs_sou_sa2, detailed = TRUE)
  l.out_bri_sa2 <- f.run_scenario(l.inputs = l.inputs_bri_sa2, detailed = TRUE)
  calculate_icers(cost = l.out_eas_sa2[["df.out"]][,"COST"], effect = l.out_eas_sa2[["df.out"]][,"QALY"], strategies = l.out_eas_sa2[["df.out"]][,"strategy"])
  calculate_icers(cost = l.out_nor_sa2[["df.out"]][,"COST"], effect = l.out_nor_sa2[["df.out"]][,"QALY"], strategies = l.out_nor_sa2[["df.out"]][,"strategy"])
  calculate_icers(cost = l.out_wes_sa2[["df.out"]][,"COST"], effect = l.out_wes_sa2[["df.out"]][,"QALY"], strategies = l.out_wes_sa2[["df.out"]][,"strategy"])
  calculate_icers(cost = l.out_sou_sa2[["df.out"]][,"COST"], effect = l.out_sou_sa2[["df.out"]][,"QALY"], strategies = l.out_sou_sa2[["df.out"]][,"strategy"])
  calculate_icers(cost = l.out_bri_sa2[["df.out"]][,"COST"], effect = l.out_bri_sa2[["df.out"]][,"QALY"], strategies = l.out_bri_sa2[["df.out"]][,"strategy"])
    # range: 78,000-124,000
  
  # sensitivity analysis 3: 10% waning
  l.inputs_eas_sa3 <- l.inputs_eas
  l.inputs_nor_sa3 <- l.inputs_nor
  l.inputs_wes_sa3 <- l.inputs_wes
  l.inputs_sou_sa3 <- l.inputs_sou
  l.inputs_bri_sa3 <- l.inputs_bri
  l.inputs_eas_sa3[["tx_waning"]] <- 0.10
  l.inputs_nor_sa3[["tx_waning"]] <- 0.10
  l.inputs_wes_sa3[["tx_waning"]] <- 0.10
  l.inputs_sou_sa3[["tx_waning"]] <- 0.10
  l.inputs_bri_sa3[["tx_waning"]] <- 0.10
  l.out_eas_sa3 <- f.run_scenario(l.inputs = l.inputs_eas_sa3, detailed = TRUE)
  l.out_nor_sa3 <- f.run_scenario(l.inputs = l.inputs_nor_sa3, detailed = TRUE)
  l.out_wes_sa3 <- f.run_scenario(l.inputs = l.inputs_wes_sa3, detailed = TRUE)
  l.out_sou_sa3 <- f.run_scenario(l.inputs = l.inputs_sou_sa3, detailed = TRUE)
  l.out_bri_sa3 <- f.run_scenario(l.inputs = l.inputs_bri_sa3, detailed = TRUE)
  calculate_icers(cost = l.out_eas_sa3[["df.out"]][,"COST"], effect = l.out_eas_sa3[["df.out"]][,"QALY"], strategies = l.out_eas_sa3[["df.out"]][,"strategy"])
  calculate_icers(cost = l.out_nor_sa3[["df.out"]][,"COST"], effect = l.out_nor_sa3[["df.out"]][,"QALY"], strategies = l.out_nor_sa3[["df.out"]][,"strategy"])
  calculate_icers(cost = l.out_wes_sa3[["df.out"]][,"COST"], effect = l.out_wes_sa3[["df.out"]][,"QALY"], strategies = l.out_wes_sa3[["df.out"]][,"strategy"])
  calculate_icers(cost = l.out_sou_sa3[["df.out"]][,"COST"], effect = l.out_sou_sa3[["df.out"]][,"QALY"], strategies = l.out_sou_sa3[["df.out"]][,"strategy"])
  calculate_icers(cost = l.out_bri_sa3[["df.out"]][,"COST"], effect = l.out_bri_sa3[["df.out"]][,"QALY"], strategies = l.out_bri_sa3[["df.out"]][,"strategy"])
    # range: 406,000-433,000
}



######################################## 5.5. ABSTRACT CTAD 2024 ########################################

if(F) {
  
  # copy inputs ICER replication
  l.inputs_ctad <- l.inputs_icer
  
  # add costs for amyloid PET testing (50% = prevalence abnormal PET [assumption]; 15% = proportion already receiving amyloid PET [assumption]; number needed to test = 1/(0.50) * (1-0.15); costs PET = 4230 [https://doi.org/10.1212/wnl.0000000000209218])
  l.inputs_ctad[["c.Tx_start"]] <- 261.10*4 + 261.10*3*0.215 + 1/(0.50)*(1-0.15)*4230 # MRI monitoring + MRI AE + amyloid PET (nnt*proportion_untested*PET_costs)
  
  # run scenario
  l.out_ctad <- f.run_scenario(l.inputs = l.inputs_ctad, detailed = TRUE)
  
  # headroom function
  f.headroom <- function(x, l.inputs, parameter) {
    l.inputs[[parameter]] <- x
    out <- f.run_scenario(l.inputs=l.inputs, detailed=FALSE)
    iNHB <- out[2,"NHB"] - out[1,"NHB"]
    return(iNHB)
  }
  # run headroom
  headroom_ctad <- optimize(f=function(x) abs(f.headroom(x, l.inputs=l.inputs_ctad, parameter="c.Tx")), interval=c(1,100000))[["minimum"]]; headroom_ctad
  # temporary for testing
  l.inputs_test <- l.inputs_ctad
  l.inputs_test$c.Tx <- 10000
  l.out_test <- f.run_scenario(l.inputs = l.inputs_test, detailed = TRUE)
  l.out_test$df.out
  calculate_icers(cost = l.out_test[["df.out"]][,"COST"], effect = l.out_test[["df.out"]][,"QALY"], strategies = l.out_test[["df.out"]][,"strategy"])
  
  # additional results
  m.result_ctad <- matrix(data = NA, nrow = ncol(l.out_ctad$l.out$soc$m.out), ncol = 3, dimnames = list(colnames(l.out_ctad$l.out$soc$m.out), c("soc","int","dif")))
  m.result_ctad[,"soc"] <- colSums(l.out_ctad$l.out$soc$m.out)
  m.result_ctad[,"int"] <- colSums(l.out_ctad$l.out$int$m.out)
  m.result_ctad[,"dif"] <- m.result_ctad[,"int"] - m.result_ctad[,"soc"]
  
  # RESULT: QALY gain
  m.result_ctad["qaly","dif"]
  
  # RESULT: range care sector costs/savings
  m.result_ctad[c("cost_hc","cost_sc"),"dif"]
  sum(m.result_ctad[c("cost_hc","cost_sc"),"dif"])
  # RESULT: diagnostic and drug sector costs
  m.result_ctad[c("cost_dx","cost_tx"),"dif"]
  sum(m.result_ctad[c("cost_dx","cost_tx"),"dif"])
  # RESULT: ICER
  calculate_icers(cost = l.out_ctad[["df.out"]][,"COST"], effect = l.out_ctad[["df.out"]][,"QALY"], strategies = l.out_ctad[["df.out"]][,"strategy"])
  
  # scenario: subcutaneous
  l.inputs_ctad_sub <- l.inputs_ctad
  l.inputs_ctad_sub[["c.Tx"]] <- 26500 + (52/2)*78.35*(1/3) # drug annual wholesale acquisition cost + treatment administration frequency * administration cost
  l.out_ctad_sub <- f.run_scenario(l.inputs = l.inputs_ctad_sub, detailed = TRUE)
  calculate_icers(cost = l.out_ctad_sub[["df.out"]][,"COST"], effect = l.out_ctad_sub[["df.out"]][,"QALY"], strategies = l.out_ctad_sub[["df.out"]][,"strategy"])
  
  # scenario: plasma amyloid (50% = prevalence abnormal PET [assumption]; 0.74 = sensitivity, 0.78 = specificity amyloid in plasma [https://doi-org.mu.idm.oclc.org/10.1212/wnl.0000000000013211]; 500 = costs plasma [https://doi-org.mu.idm.oclc.org/10.1212/wnl.0000000000209218])
  l.inputs_ctad_pla <- l.inputs_ctad
  l.inputs_ctad_pla[["c.Tx_start"]] <- 261.10*4 + 261.10*3*0.215 + 1/(0.50)*500 + (0.74*0.50 + (1-0.78)*(1-0.50))*(1-0.15)*4230 # MRI monitoring + MRI AE + plasma (nnt*plasma_cost) + amyloid PET (proportion plasma positive * proportion_untested * PET_cost)
  l.out_ctad_pla <- f.run_scenario(l.inputs = l.inputs_ctad_pla, detailed = TRUE)
  calculate_icers(cost = l.out_ctad_pla[["df.out"]][,"COST"], effect = l.out_ctad_pla[["df.out"]][,"QALY"], strategies = l.out_ctad_pla[["df.out"]][,"strategy"])
  
  # scenario: subgroup selection: starting population mci/mil
  l.inputs_ctad_mci <- l.inputs_ctad
  l.inputs_ctad_mci[["p.starting_state_mci"]] <- 1
  l.out_ctad_mci <- f.run_scenario(l.inputs = l.inputs_ctad_mci, detailed = TRUE)
  calculate_icers(cost = l.out_ctad_mci[["df.out"]][,"COST"], effect = l.out_ctad_mci[["df.out"]][,"QALY"], strategies = l.out_ctad_mci[["df.out"]][,"strategy"])
  l.inputs_ctad_mil <- l.inputs_ctad
  l.inputs_ctad_mil[["p.starting_state_mci"]] <- 0
  l.out_ctad_mil <- f.run_scenario(l.inputs = l.inputs_ctad_mil, detailed = TRUE)
  calculate_icers(cost = l.out_ctad_mil[["df.out"]][,"COST"], effect = l.out_ctad_mil[["df.out"]][,"QALY"], strategies = l.out_ctad_mil[["df.out"]][,"strategy"])
  
  # scenario: subgroup selection: ApoE4 noncarrier/carrier
  l.inputs_ctad_noncarrier <- l.inputs_ctad
  l.inputs_ctad_noncarrier[["rr.tx_mci_mil"]] <- l.inputs_ctad_noncarrier[["rr.tx_mil_mod"]] <- l.inputs_ctad_noncarrier[["rr.tx_mil_sev"]] <- 0.59
  l.out_ctad_noncarrier <- f.run_scenario(l.inputs = l.inputs_ctad_noncarrier, detailed = TRUE)
  calculate_icers(cost = l.out_ctad_noncarrier[["df.out"]][,"COST"], effect = l.out_ctad_noncarrier[["df.out"]][,"QALY"], strategies = l.out_ctad_noncarrier[["df.out"]][,"strategy"])
  l.inputs_ctad_carrier <- l.inputs_ctad
  l.inputs_ctad_carrier[["rr.tx_mci_mil"]] <- l.inputs_ctad_noncarrier[["rr.tx_mil_mod"]] <- l.inputs_ctad_noncarrier[["rr.tx_mil_sev"]] <- 0.79
  l.out_ctad_carrier <- f.run_scenario(l.inputs = l.inputs_ctad_carrier, detailed = TRUE)
  calculate_icers(cost = l.out_ctad_carrier[["df.out"]][,"COST"], effect = l.out_ctad_carrier[["df.out"]][,"QALY"], strategies = l.out_ctad_carrier[["df.out"]][,"strategy"])
  
  # scenario: assume sustained effect waning
  l.inputs_ctad_suswane <- l.inputs_ctad
  l.inputs_ctad_suswane[["rr.tx_mci_mil_dis"]] <- l.inputs_ctad_suswane[["rr.tx_mci_mil"]]
  l.inputs_ctad_suswane[["rr.tx_mil_mod_dis"]] <- l.inputs_ctad_suswane[["rr.tx_mil_mod"]]
  l.inputs_ctad_suswane[["rr.tx_mil_sev_dis"]] <- l.inputs_ctad_suswane[["rr.tx_mil_sev"]]
  l.inputs_ctad_suswane[["p.tx_discontinuation2"]] <- 1
  l.inputs_ctad_suswane[["tx_discontinuation2_begin"]] <- 2
  l.inputs_ctad_suswane[["tx_waning_dis"]] <- 0.30
  l.out_ctad_suswane <- f.run_scenario(l.inputs = l.inputs_ctad_suswane, detailed = TRUE)
  calculate_icers(cost = l.out_ctad_suswane[["df.out"]][,"COST"], effect = l.out_ctad_suswane[["df.out"]][,"QALY"], strategies = l.out_ctad_suswane[["df.out"]][,"strategy"])
  
  # scenario: assume sustained effect no waning
  l.inputs_ctad_susnowane <- l.inputs_ctad_suswane
  l.inputs_ctad_susnowane[["tx_waning_dis"]] <- 0
  l.out_ctad_susnowane <- f.run_scenario(l.inputs = l.inputs_ctad_susnowane, detailed = TRUE)
  calculate_icers(cost = l.out_ctad_susnowane[["df.out"]][,"COST"], effect = l.out_ctad_susnowane[["df.out"]][,"QALY"], strategies = l.out_ctad_susnowane[["df.out"]][,"strategy"])
  
  # aggregate
  t.names <- names(colSums(l.out_ctad$l.out$soc$m.out))
  a.combine <- array(data = NA, dim = c(29,3,7), dimnames = list(t.names,c("soc","int","dif"),c("base","sub","pla","mci","noncarrier","suswane","susnowane")))
  a.combine[,"soc","base"] <- colSums(l.out_ctad$l.out$soc$m.out)
  a.combine[,"int","base"] <- colSums(l.out_ctad$l.out$int$m.out)
  a.combine[,"soc","sub"] <- colSums(l.out_ctad_sub$l.out$soc$m.out)
  a.combine[,"int","sub"] <- colSums(l.out_ctad_sub$l.out$int$m.out)
  a.combine[,"soc","pla"] <- colSums(l.out_ctad_pla$l.out$soc$m.out)
  a.combine[,"int","pla"] <- colSums(l.out_ctad_pla$l.out$int$m.out)
  a.combine[,"soc","mci"] <- colSums(l.out_ctad_mci$l.out$soc$m.out)
  a.combine[,"int","mci"] <- colSums(l.out_ctad_mci$l.out$int$m.out)
  a.combine[,"soc","noncarrier"] <- colSums(l.out_ctad_noncarrier$l.out$soc$m.out)
  a.combine[,"int","noncarrier"] <- colSums(l.out_ctad_noncarrier$l.out$int$m.out)
  a.combine[,"soc","suswane"] <- colSums(l.out_ctad_suswane$l.out$soc$m.out)
  a.combine[,"int","suswane"] <- colSums(l.out_ctad_suswane$l.out$int$m.out)
  a.combine[,"soc","susnowane"] <- colSums(l.out_ctad_susnowane$l.out$soc$m.out)
  a.combine[,"int","susnowane"] <- colSums(l.out_ctad_susnowane$l.out$int$m.out)
  a.combine[,"dif",] <- a.combine[,"int",] - a.combine[,"soc",]
  
  # plots
  a.combine_pos <- a.combine_neg <- a.combine[c("mci","mil","mod","sev"),c("dif"),c(2,3,4,5,6,7,1)]
  a.combine_pos[a.combine_pos<0] <- 0
  a.combine_neg[a.combine_neg>=0] <- 0
  par(mar=c(5, 6, 4, 1), xpd=TRUE)
  barplot(height = a.combine_pos, col = c("green","yellow","orange","red"), density = c(NA), beside = FALSE, horiz = TRUE, main = "difference in person-years in state", xlim = c(-1,1), las=1)
  barplot(height = a.combine_neg, col = c("green","yellow","orange","red"), density = c(NA), beside = FALSE, horiz = TRUE, ylab = "", main = "", add = TRUE, axes = FALSE, axisnames = F)
  
  barplot(height = a.combine["qaly","dif",c(2,3,4,5,6,7,1)], col = "lightgreen", density = c(NA), beside = FALSE, horiz = TRUE, main = "difference in total QALYs", las=1)
  barplot(height = a.combine["cost","dif",c(2,3,4,5,6,7,1)], col = "lightblue", density = c(NA), beside = FALSE, horiz = TRUE, main = "difference in total costs", las=1)
  
}
