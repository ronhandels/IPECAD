
######################################## INFORMATION ########################################

# see readme.md for details



######################################## TECHNICAL PREPARATION ########################################

cat("\014") # clear console
rm(list = ls()) # clear environment
# install.packages("dampack") # remove '#' at beginning of the line and run once to install this package
library(dampack) # load package
setwd("~/GitHub/IPECAD")



######################################## 1. INPUTS ########################################

######################################## 1.1. ESTIMATED ########################################

# U.S. general population life table 2019 from ssa.gov
## import life table and select only men/women and drop age column (make sure age corresponds to row number, i.e., start with age = 1)
m.lifetable_US_2019 <- as.matrix(read.csv(file="life_tables/lifetable_US_2019_ssa.csv", header=TRUE))[,c("male","female")]
## convert probability to rate
m.mortality_rate_US_2019 <- -log(1-(m.lifetable_US_2019))
## weight rate 52% female 48% male
m.mortality_rate_US_2019 <- cbind(m.mortality_rate_US_2019, weighted=NA)
m.mortality_rate_US_2019[,"weighted"] <- m.mortality_rate_US_2019[,"male"] * 0.48 + m.mortality_rate_US_2019[,"female"] * 0.52

# U.S. general population life table 2016 from cdc.gov
m.lifetable_US_2016 <- as.matrix(read.csv(file="life_tables/lifetable_US_2016.csv", header=TRUE))[,c("male","female","total")]
m.mortality_rate_US_2016 <- -log(1-(m.lifetable_US_2016))
m.mortality_rate_US_2016 <- cbind(m.mortality_rate_US_2016, weighted=NA)
m.mortality_rate_US_2016[,"weighted"] <- m.mortality_rate_US_2016[,"male"] * (1-0.446) + m.mortality_rate_US_2016[,"female"] * 0.446




######################################## 1.2. MODEL INPUTS LIST ########################################

######################################## 1.2.1. INPUTS: CROSS-VALIDATION ICER ########################################

# input parameters
l.inputs_icer <- list(
  v.names_state = c("mcion_c","mciof_c","milon_c","milof_c","mod_c","sev_c","mci_i","mil_i","mod_i","sev_i","dth"), # disease states: mci = mild cognitive impairment; mil = mild dementia; mod = moderate dementia; sev = severe dementia; dth = dead; x_i = living in institutional setting (without '_i' = living in community)
  v.names_strat = c("soc","int"), 
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

# input parameters
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

# The model is run using 2 functions: run strategy and run scenario. The second (scenario) includes a loop over all strategies by calling them one by one. This is done in the following steps (annotated in the code): 
# A: function to run a scenario
# B: prepare and initialize objects to store scenario and strategy outcomes
# C: run each strategy in a loop
# D: prepare inputs to be used in each strategy
# E: run preparations specific for the intervention strategy
# F: store newly created or updated inputs to be used for the function to run a single strategy
# G: function to run the strategy
  # G1: prepare transition probability matrix
  # G2: some checks
  # G3: initialize objects to store strategy outcomes
  # G4: starting state
  # G5: markov multiplication by looping over cycles
  # G6: multiply states with utility and cost estimates
  # G7: half-cycle correction
  # G8: discount QALYs and costs
  # G9: store outcomes to be wrapped up by the 'run scenario' function
# H: store strategy results
# I: add strategy results to scenario outcomes



######################################## 2.1. RUN STRATEGY (STEP G1-G8) ########################################

# run strategy (STEP G: function for running a strategy)
f.run_strategy <- function(l.inputs) {
  with(as.list(l.inputs), {
    
    # initialize time-dependent TP matrix (STEP G1: prepare transition probability matrix)
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
    
    # check TPs sum to 1 for each cycle (STEP G2: some checks)
    for(i in v.names_state) {
      temp1 <- colSums(a.TP[i,,])
      if(!isTRUE(all.equal(current = temp1, target = rep(1,n.cycle), tolerance = 1e-10))) stop(paste("TPs for",i,"do not add up to 1"))
    }
    # !!! TO-DO: check TPs are within 0-1 range
    
    
    # initialize state trace (STEP G3: initialize objects to store strategy outcomes)
    m.trace <- matrix(data = NA, nrow = n.cycle, ncol = n.state, dimnames = list(NULL,v.names_state))
    
    # initialize output table
    m.out <- matrix(data = NA, nrow = n.cycle, ncol = 29, dimnames = list(NULL, c(colnames(m.trace),"qaly_pt","qaly_ic","cost_dx","cost_tx","cost_hc","cost_sc","cost_ic","mci","mil","mod","sev","commun","instit","ontx","ly","qaly","cost","nhb")))
    
    # set starting state distribution (STEP G4: starting state)
    m.trace[1,] <- m.trace1
    
    # markov multiplication (STEP G5: markov multiplication by looping over cycles)
    for(t in 1:(n.cycle-1)) {
      m.trace[t+1,] <- m.trace[t,] %*% a.TP[,,t]
    }
    
    # primary economic outputs (STEP G6: multiply states with utility and cost estimates)
    m.out[,colnames(m.trace)] <- m.trace
    m.out[,"qaly_pt"] <- m.trace %*% c(u.mci_pt, u.mci_pt, u.mil_pt, u.mil_pt, u.mod_pt, u.sev_pt, u.mci_pt_i, u.mil_pt_i, u.mod_pt_i, u.sev_pt_i, 0) # must match order of states
    m.out[,"qaly_ic"] <- m.trace %*% c(u.mci_ic, u.mci_ic, u.mil_ic, u.mil_ic, u.mod_ic, u.sev_ic, u.mci_ic_i, u.mil_ic_i, u.mod_ic_i, u.sev_ic_i, 0) # must match order of states
    m.out[,"cost_dx"] <- 0
    m.out[,"cost_tx"] <- m.trace %*% c(c.Tx    , 0       , c.Tx    , 0       , 0       , 0       , 0         , 0         , 0         , 0         , 0) # must match order of states
    m.out[,"cost_hc"] <- m.trace %*% c(c.mci_hc, c.mci_hc, c.mil_hc, c.mil_hc, c.mod_hc, c.sev_hc, c.mci_hc_i, c.mil_hc_i, c.mod_hc_i, c.sev_hc_i, 0) # must match order of states
    m.out[,"cost_sc"] <- m.trace %*% c(c.mci_sc, c.mci_sc, c.mil_sc, c.mil_sc, c.mod_sc, c.sev_sc, c.mci_sc_i, c.mil_sc_i, c.mod_sc_i, c.sev_sc_i, 0) # must match order of states
    m.out[,"cost_ic"] <- m.trace %*% c(c.mci_ic, c.mci_ic, c.mil_ic, c.mil_ic, c.mod_ic, c.sev_ic, c.mci_ic_i, c.mil_ic_i, c.mod_ic_i, c.sev_ic_i, 0) # must match order of states
    
    # half-cycle correction (STEP G7: apply half-cycle correction)
    if(half_cycle_correction) {
      for (j in colnames(m.out)) {
        for (i in 1:(n.cycle-1)) {
          m.out[i,j]   <- (m.out[i,j]   + m.out[i+1,j])   * 0.5
        }
      }
      m.out <- m.out[-n.cycle,]
    }
    
    # add additional outcomes at cycle 1
    if(strat=="int") {
      m.out[,"cost_dx"][1] <- m.out[,"cost_dx"][1] + c.Tx_start
      m.out[,"qaly_pt"][1] <- m.out[,"qaly_pt"][1] + u.Tx_start
    }
    
    # define vector for discounting QALYs and costs (STEP G8: apply discounting)
    n <- ifelse(test=half_cycle_correction, yes=2, no=1)
    v.discount_EFFECT <- 1 / (( 1 + (discount_EFFECT)) ^ (0 : (n.cycle-n)))
    v.discount_QALY   <- 1 / (( 1 + (discount_QALY))   ^ (0 : (n.cycle-n)))
    v.discount_COST   <- 1 / (( 1 + (discount_COST))   ^ (0 : (n.cycle-n)))
    
    # apply discounting
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
    
    # store strategy-specific output (STEP G9: store outcomes to be wrapped up by the 'run scenario' function)
    return(list(
      a.TP = a.TP, 
      m.trace = m.trace, 
      m.out = m.out
    ))
  }
  )
}



######################################## 2.2. RUN SCENARIO (STEP A-I) ########################################

# run scenario (STEP A: function for running a scenario, mainly to loop over strategies)
f.run_scenario <- function(l.inputs, detailed=FALSE) {
  with(as.list(l.inputs), { # the 'with' functions enables to call the items from the list without having to refer to the list each time (one can use 'age_start' instead l.inputs[["age_start"]])
    
    # store counters (STEP B: prepare and initialize objects to store scenario and strategy outcomes)
    n.state <- length(v.names_state) # number of states
    n.strat <- length(v.names_strat) # number of strategies

    # initialize output dataframe (create an empty dataframe to store outcomes of a scenario)
    df.out <- data.frame(
      strategy = v.names_strat,
      QALY = numeric(n.strat),
      COST = numeric(n.strat),
      LY = numeric(n.strat),
      NHB = numeric(n.strat),
      row.names = v.names_strat, 
      stringsAsFactors = FALSE
    )
    
    # initialize output list (create empty list to store outcomes of each strategy)
    l.out <- vector(mode = "list", length = 0)
    
    # loop over strategies (STEP C: run each strategy in a loop)
    for(strat in v.names_strat) {
      
      # convert time-independent transitions to vector of transitions (STEP D: prepare inputs to be used in each strategy)
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
      
      # probability of remaining in the same state (calculated as 1 minus transitions to other states, conditional on survival)
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
      
      # strategy-specific inputs (STEP E: run preparations specific for the intervention strategy)
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
      
      # list inputs for running each strategy (STEP F: store newly created or updated inputs to be used for the function to run a single strategy)
      l.inputs_strategy <- c(l.inputs, list(
        strat = strat, 
        n.state = n.state, 
        n.strat = n.strat, 
        n.cycle = n.cycle, 
        v.p.mci_mil = v.p.mci_mil, 
        v.p.mci_mod = v.p.mci_mod, 
        v.p.mci_sev = v.p.mci_sev, 
        v.p.mil_mci = v.p.mil_mci, 
        v.p.mil_mod = v.p.mil_mod, 
        v.p.mil_sev = v.p.mil_sev, 
        v.p.mod_mil = v.p.mod_mil, 
        v.p.mod_sev = v.p.mod_sev, 
        v.p.sev_mil = v.p.sev_mil, 
        v.p.sev_mod = v.p.sev_mod, 
        v.p.mci_mci = v.p.mci_mci, 
        v.p.mil_mil = v.p.mil_mil, 
        v.p.mod_mod = v.p.mod_mod, 
        v.p.sev_sev = v.p.sev_sev, 
        v.p.mcion_mci = v.p.mcion_mci, 
        v.p.mcion_mil = v.p.mcion_mil, 
        v.p.mcion_mod = v.p.mcion_mod, 
        v.p.mcion_sev = v.p.mcion_sev, 
        v.p.milon_mci = v.p.milon_mci, 
        v.p.milon_mil = v.p.milon_mil, 
        v.p.milon_mod = v.p.milon_mod,
        v.p.milon_sev = v.p.milon_sev, 
        v.p.discontinuation = v.p.discontinuation, 
        v.r.dth = v.r.dth, 
        m.trace1 = m.trace1
      ))
      
      # run strategy (STEP G: run the strategy)
      l.strat <- f.run_strategy(l.inputs_strategy)
      
      # store strategy-specific output (STEP H store strategy results)
      l.out[[strat]] <- l.strat
      
      # store output (STEP I add strategy results to scenario outcomes)
      m.out <- l.strat[["m.out"]] # extract cycle-specific outcomes
      df.out[strat,"strategy"] <- strat # store strategy name
      df.out[strat,"QALY"] <- sum(m.out[,"qaly"]) # calculate total QALYs and store them
      df.out[strat,"COST"] <- sum(m.out[,"cost"]) # calculate total costs and store them
      df.out[strat,"LY"]   <- sum(m.out[,"ly"]) # calculate total live years and store them
      df.out[strat,"NHB"]  <- sum(m.out[,"nhb"]) # calculate total QALYs and store them
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
  )
}



######################################## 2.3. PREPARE RESULTS ########################################

f.result <- function(l.out_scenario, within) {
  
  m.out_soc <- l.out_scenario$l.out$soc$m.out
  m.out_int <- l.out_scenario$l.out$int$m.out
  
  # initialize matrix for summary outcomes
  m.result <- matrix(
    data = NA, 
    nrow = ncol(m.out_soc), 
    ncol = 11, 
    dimnames = list(
      colnames(m.out_soc),
      c("soc","int","dif","dif_relative","soc_within","int_within","soc_extrapolate","int_extrapolate","dif_within","dif_extrapolate","dif_p_extrapolate")
    )
  )

  # fill matrix
  for(i in rownames(m.result)) {
    m.result[i,"soc"] <- sum(m.out_soc[,i])
    m.result[i,"int"] <- sum(m.out_int[,i])
    m.result[i,"soc_within"] <- sum(m.out_soc[1:within,i])
    m.result[i,"int_within"] <- sum(m.out_int[1:within,i])
    m.result[i,"soc_extrapolate"] <- sum(m.out_soc[(within+1):nrow(m.out_soc),i])
    m.result[i,"int_extrapolate"] <- sum(m.out_int[(within+1):nrow(m.out_soc),i])
  }
  
  # add absolute, relative and proportional difference
  round(m.result,1)
  m.result[,"dif"] <- m.result[,"int"] - m.result[,"soc"]
  m.result[,"dif_relative"] <- m.result[,"dif"] / m.result[,"soc"]
  m.result[,"dif_within"] <- m.result[,"int_within"] - m.result[,"soc_within"]
  m.result[,"dif_extrapolate"] <- m.result[,"int_extrapolate"] - m.result[,"soc_extrapolate"]
  m.result[,"dif_p_extrapolate"] <- m.result[,"dif_extrapolate"] / m.result[,"dif"]
  round(m.result, digits=2)
    # note: proportional difference is invalid for outcomes with negative values
  
  # return
  return(m.result)
  
}



######################################## 3. MODEL CALIBRATION ########################################
# n/a



######################################## 4. VALIDATION ########################################

######################################## 4.1. COMPARE TO PREVIOUS VERSION ########################################

if(F) {
  # internal validation to previous version
  l.inputs_val_int <- list(
    v.names_state = c("mcion","mciof","milon","milof","mod","sev","mci_i","mil_i","mod_i","sev_i","dth"), # disease states: mci = mild cognitive impairment; mil = mild dementia; mod = moderate dementia; sev = severe dementia; dth = dead; x_i = living in institutional setting (without '_i' = living in community)
    v.names_strat = c("soc","int"), # strategies: soc = standard of care strategy; int = intervention strategy
    age_start = 70, 
    n.cycle = 29, 
    sex = "total", 
    p.mci_mil = 0.21, 
    p.mci_mod = 0, 
    p.mci_sev = 0, 
    p.mil_mci = 0, 
    p.mil_mod = 0.293, 
    p.mil_sev = 0.001, 
    p.mod_mil = 0.087, 
    p.mod_sev = 0.109, 
    p.sev_mil = 0.000, 
    p.sev_mod = 0.196, 
    p.mci_i = 0, 
    p.mil_i = 0.038, 
    p.mod_i = 0.110, 
    p.sev_i = 0.259, 
    m.r.mortality = m.mortality_rate_US, 
    hr.mort_mci = 1, 
    hr.mort_verymilddem = 1.82, 
    hr.mort_mil = 1.318, 
    hr.mort_mod = 2.419, 
    hr.mort_sev = 4.267, 
    rr.tx_mci_mil = 0.75, 
    rr.tx_mci_mod = 1, 
    rr.tx_mci_sev = 1, 
    rr.tx_mil_mod = 0.75, 
    rr.tx_mci_mil_dis = 1, 
    rr.tx_mci_mod_dis = 1, 
    rr.tx_mci_sev_dis = 1, 
    rr.tx_mil_mod_dis = 1, 
    p.tx_discontinuation1 = 0, 
    p.tx_discontinuation2 = 0.1, 
    tx_discontinuation2_begin = 2, 
    tx_waning = 0.05, 
    tx_waning_dis = 0, 
    tx_duration = 7, 
    p.starting_state_mci = 1, 
    u.mci_pt = 0.73, 
    u.mil_pt = 0.69, 
    u.mod_pt = 0.53, 
    u.sev_pt = 0.38, 
    u.mci_ic = 0, 
    u.mil_ic = 0, 
    u.mod_ic = 0, 
    u.sev_ic = 0, 
    c.mci_hc = 1254 * 12, 
    c.mil_hc = 1471 * 12, 
    c.mod_hc = 1958 * 12, 
    c.sev_hc = 2250 * 12, 
    c.mci_hc_i = 1254 * 12, 
    c.mil_hc_i = 1471 * 12, 
    c.mod_hc_i = 1958 * 12, 
    c.sev_hc_i = 2250 * 12, 
    c.mci_sc = 222 * 12, 
    c.mil_sc = 410 * 12, 
    c.mod_sc = 653 * 12, 
    c.sev_sc = 1095 * 12, 
    c.mci_sc_i = 8762 * 12, 
    c.mil_sc_i = 8762 * 12, 
    c.mod_sc_i = 8762 * 12, 
    c.sev_sc_i = 8762 * 12, 
    c.mci_ic = 0, 
    c.mil_ic = 0, 
    c.mod_ic = 0, 
    c.sev_ic = 0, 
    c.mci_ic_i = 0, 
    c.mil_ic_i = 0, 
    c.mod_ic_i = 0, 
    c.sev_ic_i = 0, 
    c.Tx = 10000, 
    c.Tx_start = 2000, 
    discount_QALY = 0.035, 
    discount_COST = 0.035, 
    wtp = 40000, 
    half_cycle_correction = FALSE
  )
  
  # run scenario
  out_val_int <- f.run_scenario(l.inputs = l.inputs_val_int, detailed = TRUE)
  out_val_int
  
  ## incremental summary outcomes
  df.he_incr_val_int <- rbind(
    out_val_int[["df.out"]][,-1], 
    incremental = out_val_int[["df.out"]]["int",-1] - out_val_int[["df.out"]]["soc",-1]
  )
  table.he_incr_val_int <- format(df.he_incr_val_int, digits=2, scientific=FALSE, big.mark=",")
  print(table.he_incr_val_int)
  
  # comment: fully replicated outcomes compared to previous version (main branch at this moment)

}



######################################## 4.1. EXTREME SCENARIOS ########################################

if(F) {

  # identical scenarios
  l.inputs_icer_extr1 <- l.inputs_icer
  l.inputs_icer_extr1[["rr.tx_mci_mil"]] <- 1
  l.inputs_icer_extr1[["rr.tx_mil_mod"]] <- 1
  l.inputs_icer_extr1[["rr.tx_mil_sev"]] <- 1
  l.inputs_icer_extr1[["u.Tx_start"]] <- 0
  l.inputs_icer_extr1[["c.Tx"]] <- 0
  l.inputs_icer_extr1[["c.Tx_start"]] <- 0
  ## run scenario and results
  l.out_icer_extr1 <- f.run_scenario(l.inputs = l.inputs_icer_extr1, detailed = TRUE)
  m.result_icer_extr1 <- f.result(l.out_scenario = l.out_icer_extr1, within = 2)
  ## print results
  round(m.result_icer_extr1,2)
  
}


######################################## 5. ANALYSIS ########################################

######################################## 5.1. CROSS-VALIDATION: ICER ########################################

if(T) {
  
  # run scenario and results
  l.out_icer <- f.run_scenario(l.inputs = l.inputs_icer, detailed = TRUE)
  m.result_icer <- f.result(l.out_scenario = l.out_icer, within = 2)
  
  # compare to publication
  print(round(m.result_icer[c("ly","qaly","cost"),c("soc","int","dif")],2))
  icer_icer <- calculate_icers(cost = l.out_icer[["df.out"]][,"COST"], effect = l.out_icer[["df.out"]][,"QALY"], strategies = l.out_icer[["df.out"]][,"strategy"])
  print(icer_icer)
  
  # additional analysis (presentation Stockholm 2024)
  round(m.result_icer,2)
  with(as.list(l.inputs_icer), {
    
    # prepare
    c.mci <- c.mci_hc*(1-p.mci_i) + c.mci_hc_i*p.mci_i + c.mci_sc*(1-p.mci_i) + c.mci_sc_i*p.mci_i + c.mci_ic*(1-p.mci_i) + c.mci_ic_i*p.mci_i
    c.mil <- c.mil_hc*(1-p.mil_i) + c.mil_hc_i*p.mil_i + c.mil_sc*(1-p.mil_i) + c.mil_sc_i*p.mil_i + c.mil_ic*(1-p.mil_i) + c.mil_ic_i*p.mil_i
    c.mod <- c.mod_hc*(1-p.mod_i) + c.mod_hc_i*p.mod_i + c.mod_sc*(1-p.mod_i) + c.mod_sc_i*p.mod_i + c.mod_ic*(1-p.mod_i) + c.mod_ic_i*p.mod_i
    c.sev <- c.sev_hc*(1-p.sev_i) + c.sev_hc_i*p.sev_i + c.sev_sc*(1-p.sev_i) + c.sev_sc_i*p.sev_i + c.sev_ic*(1-p.sev_i) + c.sev_ic_i*p.sev_i
    u.mci <- (u.mci_pt+u.mci_ic)*(1-p.mci_i) + (u.mci_pt_i+u.mci_ic_i)*p.mci_i
    u.mil <- (u.mil_pt+u.mil_ic)*(1-p.mil_i) + (u.mil_pt_i+u.mil_ic_i)*p.mil_i
    u.mod <- (u.mod_pt+u.mod_ic)*(1-p.mod_i) + (u.mod_pt_i+u.mod_ic_i)*p.mod_i
    u.sev <- (u.sev_pt+u.sev_ic)*(1-p.sev_i) + (u.sev_pt_i+u.sev_ic_i)*p.sev_i
    
    # outcome 1
    c.Tx_start # diagnostic costs + side effects
    c.Tx_tot <- sum(m.result_icer[c("mci","mil"),"soc"]) * c.Tx * (1-p.tx_discontinuation1) # treatment cost
    u.Tx_start # treatment side effects
    c.mci_dif <- m.result_icer["mci","dif"] * c.mci
    c.mil_dif <- m.result_icer["mil","dif"] * c.mil
    c.mod_dif <- m.result_icer["mod","dif"] * c.mod
    c.sev_dif <- m.result_icer["sev","dif"] * c.sev
    u.mci_dif <- m.result_icer["mci","dif"] * u.mci
    u.mil_dif <- m.result_icer["mil","dif"] * u.mil
    u.mod_dif <- m.result_icer["mod","dif"] * u.mod
    u.sev_dif <- m.result_icer["sev","dif"] * u.sev
    c_dif1 <- c.Tx_start + c.Tx_tot + c.mci_dif + c.mil_dif + c.mod_dif + c.sev_dif
    qaly_dif1 <- u.Tx_start + u.mci_dif + u.mil_dif + u.mod_dif + u.sev_dif
    
    # outcome 2
    c.mci_dif2 <- m.result_icer["mci","dif"] * (c.mil - c.mci)
    c.mil_dif2 <- m.result_icer["mil","dif"] * (c.mod - c.mil)
    c.ly_dif <- m.result_icer["ly","dif"] * ((c.mci + c.mil)/2)
    u.mci_dif2 <- m.result_icer["mci","dif"] * (u.mil - u.mci)
    u.mil_dif2 <- m.result_icer["mil","dif"] * (u.mod - u.mil)
    u.ly <- m.result_icer["ly","dif"] * ((u.mci + u.mil)/2)
    c_dif2 <- c.Tx_start + c.Tx_tot + -c.mci_dif2 + -c.mil_dif2 + c.ly_dif
    qaly_dif2 <- u.Tx_start + -u.mci_dif2 + -u.mil_dif2 + u.ly
    
    # return
    return(list(
      c_dif1 = c_dif1, 
      qaly_dif1 = qaly_dif1, 
      icer1 = c_dif1/qaly_dif1, 
      c_dif2 = c_dif2, 
      qaly_dif2 = qaly_dif2, 
      icer2 = c_dif2/qaly_dif2
    ))
  }
  )
  
  # standard tables/plots
  
  if(F) {
    
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
    
    # table: summary results
    table.summary_result_data <- m.result_icer[c("mci","mil","mod","sev","ly","qaly","cost"),c("soc","int","dif")]
    table.summary_result <- format(
      table.summary_result_data, 
      digits=2, 
      scientific=FALSE, 
      big.mark=","
    )
    
    print(round(table.summary_result_data,2))
    print(table.summary_result)
    
    # plot: time in state
    plot.timestate_data <- m.result_icer[c("mci","mil","mod","sev"),c("int","soc")]
    
    windows(width=7, height=4, pointsize=12)
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
    text(x=c(0,cumsum(plot.timestate_data[1:3,"int"])), y=0.5, labels=c("mci","mil","mod","sev"), pos=4)
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
    
    windows(width=7, height=7, pointsize=12)
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
    
    # plot: incremental cost-effectiveness plane
    par(mar=c(5, 4, 4, 1)+0.1, xpd=FALSE)
    windows(width=7, height=7, pointsize=12)
    print(plot(icer_icer, label="all"))
    
    # table: incremental cost-effectiveness ratio
    print(as.data.frame(icer_icer))
    
    # figure: annual cost difference by sector over time
    m.cost_incr.pos <- m.cost_incr.neg <- l.out_icer[["l.out"]][["int"]][["m.out"]][,c("cost_dx","cost_tx","cost_hc","cost_sc","cost_ic")] - l.out_icer[["l.out"]][["soc"]][["m.out"]][,c("cost_dx","cost_tx","cost_hc","cost_sc","cost_ic")] # split positive and negative
    m.cost_incr.pos[m.cost_incr.pos<0] <- 0
    m.cost_incr.neg[m.cost_incr.neg>=0] <- 0
    
    windows(width=7, height=7, pointsize=12)
    par(mar=c(5, 4, 4, 2)+0.1, xpd=FALSE)
    barplot(
      height = t(m.cost_incr.pos),
      beside = F,
      xlab = "time (years)",
      ylab = "annual incremental costs",
      ylim = c(
        min(m.cost_incr.neg) + min(m.cost_incr.neg)*0.10, 
        max(m.cost_incr.pos) + max(m.cost_incr.pos)*0.10
        ),
      col = rainbow(5), 
      names.arg = 1:nrow(m.cost_incr.pos),
      main = "costs by sector over time"
    )
    barplot(
      height = t(m.cost_incr.neg),
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
  m.result_adace <- f.result(l.out_scenario = l.out_adace, within = 2)
  
  # compare to publication
  print(round(m.result_adace[c("ly","qaly","cost"),c("soc","int","dif")],2))
  icer_adace <- calculate_icers(cost = l.out_adace[["df.out"]][,"COST"], effect = l.out_adace[["df.out"]][,"QALY"], strategies = l.out_adace[["df.out"]][,"strategy"])
  icer_adace
  
}



######################################## 5.3. UNCERTAINTY SCENARIOS ########################################

if(T) {
  
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
  
  ## initialize outcomes table
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
    l.out_scenario_prep[[i]] <- f.result(l.out_scenario = l.out_scenario[[i]], within = 2)
    if(l.out_scenario_prep[[i]]["mci","dif"]<0) stop("lower person years")
    if(l.out_scenario_prep[[i]]["mil","dif"]<0) stop("lower person years")
    m.table1[i,"mcimil"] <- sum(l.out_scenario_prep[[i]][c("mci","mil"),"dif"])
    m.table1[i,"ly"] <- l.out_scenario_prep[[i]]["ly","dif"]
    m.table1[i,"qaly"] <- l.out_scenario_prep[[i]]["qaly","dif"]
    m.table1[i,"cost_dxtx"] <- sum(l.out_scenario_prep[[i]][c("cost_dx","cost_tx"),"dif"])
    m.table1[i,"cost_care"] <- sum(l.out_scenario_prep[[i]][c("cost_hc","cost_sc","cost_ic"),"dif"])
    m.table1[i,"nhb"] <- sum(l.out_scenario_prep[[i]]["nhb","dif"])
    m.table1[i,"icer"] <- calculate_icers(cost = l.out_scenario[[i]][["df.out"]][,"COST"], effect = l.out_scenario[[i]][["df.out"]][,"QALY"], strategies = l.out_scenario[[i]][["df.out"]][,"strategy"])[2,"ICER"]
  }
  print(round(m.table1,2))
  write.table(x = round(m.table1,2), file = "clipboard", sep = "\t", row.names = FALSE, col.names = FALSE)
  
  # proportion results extrapolated
  round(m.result_icer[c("mci","mil","ly","qaly","cost_dx","cost_tx","cost_hc","cost_sc","cost_ic"),c("dif_within","dif_extrapolate","dif_p_extrapolate")],2)
  
  # additional analysis (presentation Stockholm 19-01-2024)
  ## inputs
  l.inputs_icer_6b <- l.inputs_icer_6
  l.inputs_icer_6b[["rr.tx_mci_mil"]] <- 1-0.161
  l.inputs_icer_6b[["rr.tx_mil_mod"]] <- 1-0.161
  l.inputs_icer_6b[["rr.tx_mil_sev"]] <- 1-0.161
  ## run scenario and results
  l.out_icer_6b <- f.run_scenario(l.inputs = l.inputs_icer_6b, detailed = TRUE)
  m.result_icer_6b <- f.result(l.out_scenario = l.out_icer_6b, within = 4)
  ## outcomes
  sum(m.result_icer_6b[c("mci","mil"),"dif"])
  m.result_icer_6b["ly","dif"]
  m.result_icer_6b["qaly","dif"]
  sum(m.result_icer_6b[c("cost_dx","cost_tx"),"dif"])
  sum(m.result_icer_6b[c("cost_hc","cost_sc","cost_ic"),"dif"])
  sum(m.result_icer_6b["nhb","dif"])
  calculate_icers(cost = l.out_icer_6b[["df.out"]][,"COST"], effect = l.out_icer_6b[["df.out"]][,"QALY"], strategies = l.out_icer_6b[["df.out"]][,"strategy"])[2,"ICER"]
  
  
  
}


######################################## 5.3.1. CYCLE TIME ########################################

if(T) {
  
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
  l.inputs_cycle2<- list(
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
  
  # run scenario and results
  l.out_cycle2 <- f.run_scenario(l.inputs = l.inputs_cycle2, detailed = TRUE)
  m.result_cycle2 <- f.result(l.out_scenario = l.out_cycle2, within = 4)
  
  # store result in sensitivity analysis table
  m.table1[12,"mcimil"] <- sum(m.result_cycle2[c("mci","mil"),"dif"])/(1/t)
  m.table1[12,"ly"] <- m.result_cycle2["ly","dif"]/(1/t)
  m.table1[12,"qaly"] <- m.result_cycle2["qaly","dif"]
  m.table1[12,"cost_dxtx"] <- sum(m.result_cycle2[c("cost_dx","cost_tx"),"dif"])
  m.table1[12,"cost_care"] <- sum(m.result_cycle2[c("cost_hc","cost_sc","cost_ic"),"dif"])
  m.table1[12,"nhb"] <- m.result_cycle2["nhb","dif"]
  m.table1[12,"icer"] <- calculate_icers(cost = l.out_cycle2[["df.out"]][,"COST"], effect = l.out_cycle2[["df.out"]][,"QALY"], strategies = l.out_cycle2[["df.out"]][,"strategy"])[2,"ICER"]
  
  # compare state trace
  tracesoc_icer   <- l.out_icer  [["l.out"]][["soc"]][["m.trace"]] # store trace
  traceint_icer   <- l.out_icer  [["l.out"]][["int"]][["m.trace"]]
  tracesoc_cycle2 <- l.out_cycle2[["l.out"]][["soc"]][["m.trace"]]
  traceint_cycle2 <- l.out_cycle2[["l.out"]][["int"]][["m.trace"]]
  round(tracesoc_icer  [c(1:5,10,20),],2) # print trace for cycles 1:5 and cycle 10 and cycle 20 (with 1-year cycle length)
  round(tracesoc_cycle2[(c(1:5,10,20)-1)*(1/t)+1,],2) # print trace for cycles 1,13,25,37,49 and cycle 109 and cycle 229 (with 1-month cycle length)
  round(traceint_icer  [c(1:5,10,20),],2) # idem
  round(traceint_cycle2[(c(1:5,10,20)-1)*(1/t)+1,],2) # idem
  plot (rowSums(tracesoc_icer  [c(1:10)            ,c("mcion_c","mciof_c","milon_c","milof_c","mci_i","mil_i")]), type="l", col="black", lty=1) # plot trace mci and mil for soc
  lines(rowSums(traceint_icer  [c(1:10)            ,c("mcion_c","mciof_c","milon_c","milof_c","mci_i","mil_i")]), type="l", col="black", lty=2)
  #lines(rowSums(l.out_icer$l.out$soc$m.out[,c("mci","mil")]), col="red") # half-cycle corrected
  #lines(rowSums(l.out_icer$l.out$int$m.out[,c("mci","mil")]), col="red", lty=2) # half-cycle corrected
  lines(rowSums(tracesoc_cycle2[(c(1:10)-1)*(1/t)+1,c("mcion_c","mciof_c","milon_c","milof_c","mci_i","mil_i")]), type="l", col="blue" , lty=1)
  lines(rowSums(traceint_cycle2[(c(1:10)-1)*(1/t)+1,c("mcion_c","mciof_c","milon_c","milof_c","mci_i","mil_i")]), type="l", col="blue" , lty=2)
  sum(tracesoc_icer  [2:3 ,c("mcion_c","mciof_c","milon_c","milof_c","mci_i","mil_i")]) # person-years in soc in first 2 years with 1-year cycle length
  sum(traceint_icer  [2:3 ,c("mcion_c","mciof_c","milon_c","milof_c","mci_i","mil_i")]) # person-years in int in first 2 years with 1-year cycle length
  sum(tracesoc_cycle2[2:25,c("mcion_c","mciof_c","milon_c","milof_c","mci_i","mil_i")])/12 # person-years in soc in first 2 years with 1-month cycle length
  sum(traceint_cycle2[2:25,c("mcion_c","mciof_c","milon_c","milof_c","mci_i","mil_i")])/12 # person-years in int in first 2 years with 1-month cycle length
  # compare health-economic outcomes
  l.out_icer[["df.out"]]
  l.out_cycle2[["df.out"]]
  icer_icer <- calculate_icers(cost = l.out_icer[["df.out"]][,"COST"], effect = l.out_icer[["df.out"]][,"QALY"], strategies = l.out_icer[["df.out"]][,"strategy"]); print(icer_icer, digits=2)
  icer_cycle2 <- calculate_icers(cost = l.out_cycle2[["df.out"]][,"COST"], effect = l.out_cycle2[["df.out"]][,"QALY"], strategies = l.out_cycle2[["df.out"]][,"strategy"]); print(icer_cycle2, digits=2)
  
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
    
    # run scenario and results
    l.out_cal <- f.run_scenario(l.inputs = l.inputs_cal, detailed=TRUE)
    m.result_cal <- f.result(l.out_scenario = l.out_cal, within = 2)
    
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
  
  # run scenario and results at treatment rr = 0.56 (applied to icer base case with cycle length 1 year)
  l.inputs_icer3 <- l.inputs_icer
  l.inputs_icer3[["rr.tx_mci_mil"]] <- 0.56
  l.inputs_icer3[["rr.tx_mil_mod"]] <- 0.56
  l.inputs_icer3[["rr.tx_mil_sev"]] <- 0.56
  l.out_icer3 <- f.run_scenario(l.inputs = l.inputs_icer3, detailed = TRUE)
  m.result_icer3 <- f.result(l.out_scenario = l.out_icer3, within = 18)
  
  # store result in sensitivity analysis table
  m.table1[3,"mcimil"] <- sum(m.result_icer3[c("mci","mil"),"dif"])
  m.table1[3,"ly"] <- m.result_icer3["ly","dif"]
  m.table1[3,"qaly"] <- m.result_icer3["qaly","dif"]
  m.table1[3,"cost_dxtx"] <- sum(m.result_icer3[c("cost_dx","cost_tx"),"dif"])
  m.table1[3,"cost_care"] <- sum(m.result_icer3[c("cost_hc","cost_sc","cost_ic"),"dif"])
  m.table1[3,"nhb"] <- m.result_icer3["nhb","dif"]
  m.table1[3,"icer"] <- calculate_icers(cost = l.out_icer3[["df.out"]][,"COST"], effect = l.out_icer3[["df.out"]][,"QALY"], strategies = l.out_icer3[["df.out"]][,"strategy"])[2,"ICER"]
  
  # copy sensitivity analysis results table to clipboard
  write.table(x = round(m.table1,2), file = "clipboard", sep = "\t", row.names = FALSE, col.names = FALSE)
  
}


######################################## 5.3.3. REPLICATION: HERRING ########################################


if(F) {
  
  # U.S. general population life table 2017
  m.lifetable_US_2017 <- as.matrix(read.csv(file="life_tables/lifetable_US_2017.csv", header=TRUE))[,c("male","female","total")]
  m.mortality_rate_US_2017 <- -log(1-(m.lifetable_US_2017))
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
  
  # run scenario and results
  l.out_herring <- f.run_scenario(l.inputs = l.inputs_herring, detailed = TRUE)
  m.result_herring <- f.result(l.out_scenario = l.out_herring, within = 2)
  
  # compare to publication
  print(round(m.result_herring[c("ly","ontx","mci"),c("soc","int","dif")],2))
  # dem
  print(round(sum(m.result_herring[c("mil","mod","sev"),"soc"]),2))
  print(round(sum(m.result_herring[c("mil","mod","sev"),"int"]),2))
  print(round(sum(m.result_herring[c("mil","mod","sev"),"int"]) - sum(m.result_herring[c("mil","mod","sev"),"soc"]),2))
  
}
