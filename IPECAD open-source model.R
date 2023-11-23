
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

# U.S. general population life table
m.lifetable_US <- as.matrix(read.csv(file="life_tables/lifetable_US.csv", header=TRUE))[,c("male","female")] # import life table and select only men/women and drop age column (make sure age corresponds to row number, i.e., start with age = 1)
m.mortality_rate_US <- -log(1-(m.lifetable_US)) # convert probability to rate



######################################## 1.2. MODEL INPUTS LIST ########################################

# input parameters per cycle (unless stated otherwise), see readme.md for a short explanation of each input parameter
l.inputs <- list(
  v.names_state = c("mcion","mciof","milon","milof","mod","sev","mci_i","mil_i","mod_i","sev_i","dth"), # disease states: mci = mild cognitive impairment; mil = mild dementia; mod = moderate dementia; sev = severe dementia; dth = dead; x_i = living in institutional setting (without '_i' = living in community)
  v.names_strat = c("soc","int"), # strategies: soc = standard of care strategy; int = intervention strategy
  age_start = 70, 
  n.cycle = 29, 
  sex = "female", 
  p.mci_mil = 0.248, 
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
  rr.tx_mci_mil_dis = 0.75, 
  rr.tx_mci_mod_dis = 1, 
  rr.tx_mci_sev_dis = 1, 
  rr.tx_mil_mod_dis = 0.75, 
  p.discontinuation1 = 0.1, 
  p.discontinuation2 = 0.5, 
  discontinuation2_start = 2, 
  tx_waning = 0.05, 
  tx_waning_dis = 0.15, 
  tx_duration = 3, 
  p.starting_state_mci = 1, 
  u.mci_pt = 0.73, 
  u.mil_pt = 0.69, 
  u.mod_pt = 0.53, 
  u.sev_pt = 0.38, 
  u.mci_ic = 0, 
  u.mil_ic = -0.036, 
  u.mod_ic = -0.07, 
  u.sev_ic = -0.086, 
  c.mci_hc = 1254 * 12, 
  c.mil_hc = 1471 * 12, 
  c.mod_hc = 1958 * 12, 
  c.sev_hc = 2250 * 12, 
  c.mci_sc = 222 * 12, 
  c.mil_sc = 410 * 12, 
  c.mod_sc = 653 * 12, 
  c.sev_sc = 1095 * 12, 
  c.mci_ic =  988 * 12, 
  c.mil_ic = 2184 * 12, 
  c.mod_ic = 3227 * 12, 
  c.sev_ic = 5402 * 12, 
  c.mci_i_hc = 1254 * 12, 
  c.mil_i_hc = 1471 * 12, 
  c.mod_i_hc = 1958 * 12, 
  c.sev_i_hc = 2250 * 12, 
  c.mci_i_sc = 8762 * 12, 
  c.mil_i_sc = 8762 * 12, 
  c.mod_i_sc = 8762 * 12, 
  c.sev_i_sc = 8762 * 12, 
  c.mci_i_ic =  435 * 12, 
  c.mil_i_ic =  961 * 12, 
  c.mod_i_ic = 1420 * 12, 
  c.sev_i_ic = 2377 * 12, 
  c.Tx = 5000, 
  c.Tx_diagnostics1 = 2000, 
  discount_EFFECT = 0, 
  discount_QALY = 0.035, 
  discount_COST = 0.035, 
  wtp = 40000, 
  half_cycle_correction = FALSE
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
    a.TP["mcion","dth",] <- 1-exp(-(v.r.dth * hr.mort_mci))
    a.TP["mciof","dth",] <- 1-exp(-(v.r.dth * hr.mort_mci))
    a.TP["milon","dth",] <- 1-exp(-(v.r.dth * hr.mort_mil * hr.mort_verymilddem))
    a.TP["milof","dth",] <- 1-exp(-(v.r.dth * hr.mort_mil * hr.mort_verymilddem))
    a.TP["mod",  "dth",] <- 1-exp(-(v.r.dth * hr.mort_mod * hr.mort_verymilddem))
    a.TP["sev",  "dth",] <- 1-exp(-(v.r.dth * hr.mort_sev * hr.mort_verymilddem))
    a.TP["dth",  "dth",] <- 1
    a.TP["mci_i","dth",] <- 1-exp(-(v.r.dth * hr.mort_mci))
    a.TP["mil_i","dth",] <- 1-exp(-(v.r.dth * hr.mort_mil * hr.mort_verymilddem))
    a.TP["mod_i","dth",] <- 1-exp(-(v.r.dth * hr.mort_mod * hr.mort_verymilddem))
    a.TP["sev_i","dth",] <- 1-exp(-(v.r.dth * hr.mort_sev * hr.mort_verymilddem))
    
    # TP matrix state: from mci-on community-setting
    a.TP["mcion","mcion",] <- v.p.mcion_mci * (1-p.mci_i) * (1-v.p.discontinuation) * (1-a.TP["mcion","dth",])
    a.TP["mcion","mciof",] <- v.p.mcion_mci * (1-p.mci_i) *    v.p.discontinuation  * (1-a.TP["mcion","dth",])
    a.TP["mcion","milon",] <- v.p.mcion_mil * (1-p.mci_i) * (1-v.p.discontinuation) * (1-a.TP["mcion","dth",])
    a.TP["mcion","milof",] <- v.p.mcion_mil * (1-p.mci_i) *    v.p.discontinuation  * (1-a.TP["mcion","dth",])
    a.TP["mcion","mod",]   <- v.p.mcion_mod * (1-p.mci_i)                           * (1-a.TP["mcion","dth",])
    a.TP["mcion","sev",]   <- v.p.mcion_sev * (1-p.mci_i)                           * (1-a.TP["mcion","dth",])
    a.TP["mcion","mci_i",] <- v.p.mcion_mci *    p.mci_i                            * (1-a.TP["mcion","dth",])
    a.TP["mcion","mil_i",] <- v.p.mcion_mil *    p.mci_i                            * (1-a.TP["mcion","dth",])
    a.TP["mcion","mod_i",] <- v.p.mcion_mod *    p.mci_i                            * (1-a.TP["mcion","dth",])
    a.TP["mcion","sev_i",] <- v.p.mcion_sev *    p.mci_i                            * (1-a.TP["mcion","dth",])

    # TP matrix state: from mci-off community-setting
    a.TP["mciof","mciof",] <- v.p.mci_mci   * (1-p.mci_i)                           * (1-a.TP["mciof","dth",])
    a.TP["mciof","milof",] <- v.p.mci_mil   * (1-p.mci_i)                           * (1-a.TP["mciof","dth",])
    a.TP["mciof","mod",]   <- v.p.mci_mod   * (1-p.mci_i)                           * (1-a.TP["mciof","dth",])
    a.TP["mciof","sev",]   <- v.p.mci_sev   * (1-p.mci_i)                           * (1-a.TP["mciof","dth",])
    a.TP["mciof","mci_i",] <- v.p.mci_mci   *    p.mci_i                            * (1-a.TP["mciof","dth",])
    a.TP["mciof","mil_i",] <- v.p.mci_mil   *    p.mci_i                            * (1-a.TP["mciof","dth",])
    a.TP["mciof","mod_i",] <- v.p.mci_mod   *    p.mci_i                            * (1-a.TP["mciof","dth",])
    a.TP["mciof","sev_i",] <- v.p.mci_sev   *    p.mci_i                            * (1-a.TP["mciof","dth",])
    
    # TP matrix state: from mild-on community-setting
    a.TP["milon","mcion",] <- v.p.milon_mci * (1-p.mil_i) * (1-v.p.discontinuation) * (1-a.TP["milon","dth",])
    a.TP["milon","mciof",] <- v.p.milon_mci * (1-p.mil_i) *    v.p.discontinuation  * (1-a.TP["milon","dth",])
    a.TP["milon","milon",] <- v.p.milon_mil * (1-p.mil_i) * (1-v.p.discontinuation) * (1-a.TP["milon","dth",])
    a.TP["milon","milof",] <- v.p.milon_mil * (1-p.mil_i) *    v.p.discontinuation  * (1-a.TP["milon","dth",])
    a.TP["milon","mod",]   <- v.p.milon_mod * (1-p.mil_i)                           * (1-a.TP["milon","dth",])
    a.TP["milon","sev",]   <- v.p.mil_sev   * (1-p.mil_i)                           * (1-a.TP["milon","dth",])
    a.TP["milon","mci_i",] <- v.p.milon_mci *    p.mil_i                            * (1-a.TP["milon","dth",])
    a.TP["milon","mil_i",] <- v.p.milon_mil *    p.mil_i                            * (1-a.TP["milon","dth",])
    a.TP["milon","mod_i",] <- v.p.milon_mod *    p.mil_i                            * (1-a.TP["milon","dth",])
    a.TP["milon","sev_i",] <- v.p.mil_sev   *    p.mil_i                            * (1-a.TP["milon","dth",])
    
    # TP matrix state: from mild-off community-setting
    a.TP["milof","mciof",] <- v.p.mil_mci   * (1-p.mil_i)                           * (1-a.TP["milof","dth",])
    a.TP["milof","milof",] <- v.p.mil_mil   * (1-p.mil_i)                           * (1-a.TP["milof","dth",])
    a.TP["milof","mod",]   <- v.p.mil_mod   * (1-p.mil_i)                           * (1-a.TP["milof","dth",])
    a.TP["milof","sev",]   <- v.p.mil_sev   * (1-p.mil_i)                           * (1-a.TP["milof","dth",])
    a.TP["milof","mci_i",] <- v.p.mil_mci   *    p.mil_i                            * (1-a.TP["milof","dth",])
    a.TP["milof","mil_i",] <- v.p.mil_mil   *    p.mil_i                            * (1-a.TP["milof","dth",])
    a.TP["milof","mod_i",] <- v.p.mil_mod   *    p.mil_i                            * (1-a.TP["milof","dth",])
    a.TP["milof","sev_i",] <- v.p.mil_sev   *    p.mil_i                            * (1-a.TP["milof","dth",])
    
    # TP matrix state: from moderate community-setting
    a.TP["mod","milof",] <- v.p.mod_mil     * (1-p.mod_i)                           * (1-a.TP["mod","dth",])
    a.TP["mod","mod",]   <- v.p.mod_mod     * (1-p.mod_i)                           * (1-a.TP["mod","dth",])
    a.TP["mod","sev",]   <- v.p.mod_sev     * (1-p.mod_i)                           * (1-a.TP["mod","dth",])
    a.TP["mod","mil_i",] <- v.p.mod_mil     *    p.mod_i                            * (1-a.TP["mod","dth",])
    a.TP["mod","mod_i",] <- v.p.mod_mod     *    p.mod_i                            * (1-a.TP["mod","dth",])
    a.TP["mod","sev_i",] <- v.p.mod_sev     *    p.mod_i                            * (1-a.TP["mod","dth",])
    
    # TP matrix state: from severe community-setting
    a.TP["sev","milof",] <- v.p.sev_mil     * (1-p.sev_i)                           * (1-a.TP["sev","dth",])
    a.TP["sev","mod",]   <- v.p.sev_mod     * (1-p.sev_i)                           * (1-a.TP["sev","dth",])
    a.TP["sev","sev",]   <- v.p.sev_sev     * (1-p.sev_i)                           * (1-a.TP["sev","dth",])
    a.TP["sev","mil_i",] <- v.p.sev_mil     *    p.sev_i                            * (1-a.TP["sev","dth",])
    a.TP["sev","mod_i",] <- v.p.sev_mod     *    p.sev_i                            * (1-a.TP["sev","dth",])
    a.TP["sev","sev_i",] <- v.p.sev_sev     *    p.sev_i                            * (1-a.TP["sev","dth",])

    # TP matrix state: from mci institutionalized-setting
    a.TP["mci_i","mci_i",] <- v.p.mci_mci                                           * (1-a.TP["mci_i","dth",])
    a.TP["mci_i","mil_i",] <- v.p.mci_mil                                           * (1-a.TP["mci_i","dth",])
    a.TP["mci_i","mod_i",] <- v.p.mci_mod                                           * (1-a.TP["mci_i","dth",])
    a.TP["mci_i","sev_i",] <- v.p.mci_sev                                           * (1-a.TP["mci_i","dth",])
  
    # TP matrix state: from mild institutionalized-setting
    a.TP["mil_i","mci_i",] <- v.p.mil_mci                                           * (1-a.TP["mil_i","dth",])
    a.TP["mil_i","mil_i",] <- v.p.mil_mil                                           * (1-a.TP["mil_i","dth",])
    a.TP["mil_i","mod_i",] <- v.p.mil_mod                                           * (1-a.TP["mil_i","dth",])
    a.TP["mil_i","sev_i",] <- v.p.mil_sev                                           * (1-a.TP["mil_i","dth",])
    
    # TP matrix state: from moderate institutionalized-setting
    a.TP["mod_i","mil_i",] <- v.p.mod_mil                                           * (1-a.TP["mod_i","dth",])
    a.TP["mod_i","mod_i",] <- v.p.mod_mod                                           * (1-a.TP["mod_i","dth",])
    a.TP["mod_i","sev_i",] <- v.p.mod_sev                                           * (1-a.TP["mod_i","dth",])
    
    # TP matrix state: from severe institutionalized-setting
    a.TP["sev_i","mil_i",] <- v.p.sev_mil                                           * (1-a.TP["sev_i","dth",])
    a.TP["sev_i","mod_i",] <- v.p.sev_mod                                           * (1-a.TP["sev_i","dth",])
    a.TP["sev_i","sev_i",] <- v.p.sev_sev                                           * (1-a.TP["sev_i","dth",])  
    
    # check TPs sum to 1 for each cycle (STEP G2: some checks)
    for(i in v.names_state) {
      temp1 <- colSums(a.TP[i,,])
      if(!isTRUE(all.equal(current = temp1, target = rep(1,n.cycle), tolerance = 1e-10))) stop(paste("TPs for",i,"do not add up to 1"))
    }
    # !!! TO-DO: check TPs are within 0-1 range
    
    
    # initialize state trace (STEP G3: initialize objects to store strategy outcomes)
    m.trace <- matrix(data = NA, nrow = n.cycle, ncol = n.state, dimnames = list(NULL,v.names_state))
    
    # initialize output table
    m.out <- matrix(data = NA, nrow = n.cycle, ncol = 22, dimnames = list(NULL, c(colnames(m.trace),"qaly_pt","qaly_ic","qaly","cost_dx","cost_tx","cost_hc","cost_sc","cost_ic","cost","ly","nhb")))
    
    # set starting state distribution (STEP G4: starting state)
    m.trace[1,] <- m.trace1
    
    # markov multiplication (STEP G5: markov multiplication by looping over cycles)
    for(t in 1:(n.cycle-1)) {
      m.trace[t+1,] <- m.trace[t,] %*% a.TP[,,t]
    }
    
    # primary economic outputs (STEP G6: multiply states with utility and cost estimates)
    m.out[,colnames(m.trace)] <- m.trace
    m.out[,"ly"]      <- m.trace %*% c(1       , 1       , 1       , 1       , 1       , 1       , 1         , 1         , 1         , 1         , 0) # must match order of states
    m.out[,"qaly_pt"] <- m.trace %*% c(u.mci_pt, u.mci_pt, u.mil_pt, u.mil_pt, u.mod_pt, u.sev_pt, u.mci_pt  , u.mil_pt  , u.mod_pt  , u.sev_pt  , 0) # must match order of states
    m.out[,"qaly_ic"] <- m.trace %*% c(u.mci_ic, u.mci_ic, u.mil_ic, u.mil_ic, u.mod_ic, u.sev_ic, u.mci_ic  , u.mil_ic  , u.mod_ic  , u.sev_ic  , 0) # must match order of states
    m.out[,"cost_dx"] <- 0
    m.out[,"cost_tx"] <- m.trace %*% c(c.Tx    , 0       , c.Tx    , 0       , 0       , 0       , 0         , 0         , 0         , 0         , 0) # must match order of states
    m.out[,"cost_hc"] <- m.trace %*% c(c.mci_hc, c.mci_hc, c.mil_hc, c.mil_hc, c.mod_hc, c.sev_hc, c.mci_i_hc, c.mil_i_hc, c.mod_i_hc, c.sev_i_hc, 0) # must match order of states
    m.out[,"cost_sc"] <- m.trace %*% c(c.mci_sc, c.mci_sc, c.mil_sc, c.mil_sc, c.mod_sc, c.sev_sc, c.mci_i_sc, c.mil_i_sc, c.mod_i_sc, c.sev_i_sc, 0) # must match order of states
    m.out[,"cost_ic"] <- m.trace %*% c(c.mci_ic, c.mci_ic, c.mil_ic, c.mil_ic, c.mod_ic, c.sev_ic, c.mci_i_ic, c.mil_i_ic, c.mod_i_ic, c.sev_i_ic, 0) # must match order of states
    
    # half-cycle correction (STEP G7: apply half-cycle correction)
    if(half_cycle_correction) {
      for (j in colnames(m.out)) {
        for (i in 1:(n.cycle-1)) {
          m.out[i,j]   <- (m.out[i,j]   + m.out[i+1,j])   * 0.5
        }
      }
    }
    m.out[n.cycle,] <- 0
    
    # add additional costs outside half-cycle correction: diagnostics
    if(strat=="int") {
      m.out[,"cost_dx"][1] <- m.out[,"cost_dx"][1] + c.Tx_diagnostics1
    }
    
    # define vector for discounting QALYs and costs (STEP G8: apply discounting)
    v.discount_EFFECT <- 1 / (( 1 + (discount_EFFECT)) ^ (0 : (n.cycle-1)))
    v.discount_QALY <- 1 / (( 1 + (discount_QALY)) ^ (0 : (n.cycle-1)))
    v.discount_COST <- 1 / (( 1 + (discount_COST)) ^ (0 : (n.cycle-1)))
    
    # apply discounting
    for(i in colnames(m.trace)) {
      m.out[,i] <- m.out[,i]*v.discount_EFFECT
    }
    for(i in c("qaly_pt","qaly_ic")) {
      m.out[,i] <- m.out[,i]*v.discount_QALY
    }
    for(i in c("cost_dx","cost_tx","cost_hc","cost_sc","cost_ic")) {
      m.out[,i] <- m.out[,i]*v.discount_COST
    }
    
    # totals
    m.out[,"qaly"] <- rowSums(m.out[,c("qaly_pt","qaly_ic")])
    m.out[,"cost"] <- rowSums(m.out[,c("cost_dx","cost_tx","cost_hc","cost_sc","cost_ic")])
    
    # calculate net health benefit
    m.out[,"nhb"] <- m.out[,"qaly"] - (m.out[,"cost"] / wtp)
    
    # store strategy-specific output (STEP G9: store outcomes to be wrapped up by the 'run scenario' function)
    return(list(
      a.TP=a.TP, 
      m.trace=m.trace, 
      m.out=m.out
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

    # initialize output matrix (create an empty dataframe to store outcomes of a scenario)
    df.out_sum <- data.frame(
      strategy = v.names_strat,
      QALY = numeric(n.strat),
      COST = numeric(n.strat),
      LY = numeric(n.strat),
      NHB = numeric(n.strat),
      row.names = v.names_strat, 
      stringsAsFactors = FALSE
    )
    
    # initialize output list (create empty list to store outcomes of each strategy)
    l.out_strategy <- vector(mode = "list", length = 0)
    
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
      v.p.mcion_mil <- v.p.mci_mil
      v.p.mcion_mod <- v.p.mci_mod
      v.p.mcion_sev <- v.p.mci_sev
      v.p.milon_mci <- v.p.mil_mci
      v.p.milon_mod <- v.p.mil_mod
      v.p.mcion_mci <- v.p.mci_mci
      v.p.milon_mil <- v.p.mil_mil
      
      # discontinuation
      v.p.discontinuation <- rep(x=0, times=n.cycle) # initialize vector
      v.p.discontinuation[1:(discontinuation2_start-1)] <- p.discontinuation1 # discontinuation up to discontinuation2 start cycle
      v.p.discontinuation[discontinuation2_start:n.cycle] <- p.discontinuation2 # discontinuation at discontinuation2 start cycle onward
      v.p.discontinuation[tx_duration:n.cycle] <- 1 # maximum treatment duration implemented as discontinuation
      
      # death (subset mortality table to obtain age- and sex-specific mortality)
      v.r.dth <- m.r.mortality[age_start:(age_start+n.cycle-1), sex]
      
      # starting states
      m.trace1 <- matrix(data=0, nrow=1, ncol=n.state, dimnames=list(NULL,v.names_state))
      m.trace1[,"mciof"] <- p.starting_state_mci
      m.trace1[,"milof"] <- 1-p.starting_state_mci
      
      # strategy-specific inputs (STEP E: run preparations specific for the intervention strategy)
      if(strat=="int") {
        
        # waning
        temp.waning <- (1-tx_waning)^(0:(n.cycle-1))
        temp.rr.tx_mci_mil <- rr.tx_mci_mil^temp.waning
        temp.rr.tx_mci_mod <- rr.tx_mci_mod^temp.waning
        temp.rr.tx_mci_sev <- rr.tx_mci_sev^temp.waning
        temp.rr.tx_mil_mod <- rr.tx_mil_mod^temp.waning
        temp.waning_dis <- (1-tx_waning_dis)^(0:(n.cycle-1))
        temp.rr.tx_mci_mil_dis <- rr.tx_mci_mil_dis^temp.waning_dis
        temp.rr.tx_mci_mod_dis <- rr.tx_mci_mod_dis^temp.waning_dis
        temp.rr.tx_mci_sev_dis <- rr.tx_mci_sev_dis^temp.waning_dis
        temp.rr.tx_mil_mod_dis <- rr.tx_mil_mod_dis^temp.waning_dis
        
        # update transition probabilities treatment effect: during treatment
        v.p.mcion_mil <- 1-exp(-(-log(1-p.mci_mil) * temp.rr.tx_mci_mil)) # convert probability to rate, then multiply with treatment relative risk, then convert to probability
        v.p.mcion_mod <- 1-exp(-(-log(1-p.mci_mod) * temp.rr.tx_mci_mod)) # idem
        v.p.mcion_sev <- 1-exp(-(-log(1-p.mci_sev) * temp.rr.tx_mci_sev)) # idem
        v.p.milon_mod <- 1-exp(-(-log(1-p.mil_mod) * temp.rr.tx_mil_mod)) # idem
        # update transition probabilities treatment effect: after discontinuation
        v.p.mci_mil <- 1-exp(-(-log(1-p.mci_mil) * temp.rr.tx_mci_mil_dis)) # convert probability to rate, then multiply with treatment relative risk, then convert to probability
        v.p.mci_mod <- 1-exp(-(-log(1-p.mci_mod) * temp.rr.tx_mci_mod_dis)) # idem
        v.p.mci_sev <- 1-exp(-(-log(1-p.mci_sev) * temp.rr.tx_mci_sev_dis)) # idem
        v.p.mil_mod <- 1-exp(-(-log(1-p.mil_mod) * temp.rr.tx_mil_mod_dis)) # idem
        
        # update transition probabilities of remaining in the same state
        v.p.mcion_mci <- 1 - v.p.mcion_mil - v.p.mcion_mod - v.p.mcion_sev
        v.p.milon_mil <- 1 - v.p.milon_mci - v.p.milon_mod - v.p.mil_sev
        v.p.mci_mci <- 1 - v.p.mci_mil - v.p.mci_mod - v.p.mci_sev
        v.p.mil_mil <- 1 - v.p.mil_mci - v.p.mil_mod - v.p.mil_sev
        
        # starting states
        m.trace1 <- matrix(data=0, nrow=1, ncol=n.state, dimnames=list(NULL,v.names_state))
        m.trace1[,"mcion"] <- p.starting_state_mci
        m.trace1[,"milon"] <- 1-p.starting_state_mci
        
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
        v.p.mcion_mil = v.p.mcion_mil, 
        v.p.mcion_mod = v.p.mcion_mod, 
        v.p.mcion_sev = v.p.mcion_sev, 
        v.p.milon_mci = v.p.milon_mci, 
        v.p.milon_mod = v.p.milon_mod, 
        v.p.mcion_mci = v.p.mcion_mci, 
        v.p.milon_mil = v.p.milon_mil, 
        v.p.discontinuation = v.p.discontinuation, 
        v.r.dth = v.r.dth, 
        m.trace1 = m.trace1
      ))
      
      # run strategy (STEP G: run the strategy)
      l.strat <- f.run_strategy(l.inputs_strategy)
      
      # store strategy-specific output (STEP H store strategy results)
      l.out_strategy[[strat]] <- l.strat
      
      # store output (STEP I add strategy results to scenario outcomes)
      m.out <- l.strat[["m.out"]] # extract cycle-specific outcomes
      df.out_sum[strat,"strategy"] <- strat # store strategy name
      df.out_sum[strat,"QALY"] <- sum(m.out[,"qaly"]) # calculate total QALYs and store them
      df.out_sum[strat,"COST"] <- sum(m.out[,"cost"]) # calculate total costs and store them
      df.out_sum[strat,"LY"]   <- sum(m.out[,"ly"]) # calculate total live years and store them
      df.out_sum[strat,"NHB"]  <- sum(m.out[,"nhb"]) # calculate total QALYs and store them
    }
    
    # return basic result
    if(!detailed) return(df.out_sum)
    
    # return detailed result
    if(detailed) {
      return(list(
        df.out_sum = df.out_sum,
        l.out_strategy = l.out_strategy
      ))
    }
  }
  )
}



######################################## 2.3. PREPARE OUTCOMES ########################################

f.prepare_outcomes <- function(l.out_scenario, n.cycles) {

  # trace with aggregated results
  m.trace_soc <- cbind(
    mci     = rowSums(l.out_scenario[["l.out_strategy"]][["soc"]][["m.out"]][,c("mcion","mciof","mci_i")]), 
    mil     = rowSums(l.out_scenario[["l.out_strategy"]][["soc"]][["m.out"]][,c("milon","milof","mil_i")]), 
    mod     = rowSums(l.out_scenario[["l.out_strategy"]][["soc"]][["m.out"]][,c("mod","mod_i")]), 
    sev     = rowSums(l.out_scenario[["l.out_strategy"]][["soc"]][["m.out"]][,c("sev","sev_i")]), 
    alv     = rowSums(l.out_scenario[["l.out_strategy"]][["soc"]][["m.out"]][,c("mcion","mciof","milon","milof","mod","sev","mci_i","mil_i","mod_i","sev_i")]), 
    dth     =         l.out_scenario[["l.out_strategy"]][["soc"]][["m.out"]][,c("dth")], 
    home    = rowSums(l.out_scenario[["l.out_strategy"]][["soc"]][["m.out"]][,c("mcion","mciof","milon","milof","mod","sev")]), 
    instit  = rowSums(l.out_scenario[["l.out_strategy"]][["soc"]][["m.out"]][,c("mci_i","mil_i","mod_i","sev_i")]), 
    ontx    = rowSums(l.out_scenario[["l.out_strategy"]][["soc"]][["m.out"]][,c("mcion","milon")]), 
    qaly_pt =         l.out_scenario[["l.out_strategy"]][["soc"]][["m.out"]][,c("qaly_pt")], 
    qaly_ic =         l.out_scenario[["l.out_strategy"]][["soc"]][["m.out"]][,c("qaly_ic")], 
    qaly    =         l.out_scenario[["l.out_strategy"]][["soc"]][["m.out"]][,c("qaly")], 
    cost_dx =         l.out_scenario[["l.out_strategy"]][["soc"]][["m.out"]][,c("cost_dx")], 
    cost_tx =         l.out_scenario[["l.out_strategy"]][["soc"]][["m.out"]][,c("cost_tx")], 
    cost_hc =         l.out_scenario[["l.out_strategy"]][["soc"]][["m.out"]][,c("cost_hc")], 
    cost_sc =         l.out_scenario[["l.out_strategy"]][["soc"]][["m.out"]][,c("cost_sc")], 
    cost_ic =         l.out_scenario[["l.out_strategy"]][["soc"]][["m.out"]][,c("cost_ic")], 
    cost    =         l.out_scenario[["l.out_strategy"]][["soc"]][["m.out"]][,c("cost")], 
    nhb     =         l.out_scenario[["l.out_strategy"]][["soc"]][["m.out"]][,c("nhb")]
  )
  m.trace_int <- cbind(
    mci     = rowSums(l.out_scenario[["l.out_strategy"]][["int"]][["m.out"]][,c("mcion","mciof","mci_i")]), 
    mil     = rowSums(l.out_scenario[["l.out_strategy"]][["int"]][["m.out"]][,c("milon","milof","mil_i")]), 
    mod     = rowSums(l.out_scenario[["l.out_strategy"]][["int"]][["m.out"]][,c("mod","mod_i")]), 
    sev     = rowSums(l.out_scenario[["l.out_strategy"]][["int"]][["m.out"]][,c("sev","sev_i")]), 
    alv     = rowSums(l.out_scenario[["l.out_strategy"]][["int"]][["m.out"]][,c("mcion","mciof","milon","milof","mod","sev","mci_i","mil_i","mod_i","sev_i")]), 
    dth     =         l.out_scenario[["l.out_strategy"]][["int"]][["m.out"]][,c("dth")], 
    home    = rowSums(l.out_scenario[["l.out_strategy"]][["int"]][["m.out"]][,c("mcion","mciof","milon","milof","mod","sev")]), 
    instit  = rowSums(l.out_scenario[["l.out_strategy"]][["int"]][["m.out"]][,c("mci_i","mil_i","mod_i","sev_i")]), 
    ontx    = rowSums(l.out_scenario[["l.out_strategy"]][["int"]][["m.out"]][,c("mcion","milon")]), 
    qaly_pt =         l.out_scenario[["l.out_strategy"]][["int"]][["m.out"]][,c("qaly_pt")], 
    qaly_ic =         l.out_scenario[["l.out_strategy"]][["int"]][["m.out"]][,c("qaly_ic")], 
    qaly    =         l.out_scenario[["l.out_strategy"]][["int"]][["m.out"]][,c("qaly")], 
    cost_dx =         l.out_scenario[["l.out_strategy"]][["int"]][["m.out"]][,c("cost_dx")], 
    cost_tx =         l.out_scenario[["l.out_strategy"]][["int"]][["m.out"]][,c("cost_tx")], 
    cost_hc =         l.out_scenario[["l.out_strategy"]][["int"]][["m.out"]][,c("cost_hc")], 
    cost_sc =         l.out_scenario[["l.out_strategy"]][["int"]][["m.out"]][,c("cost_sc")], 
    cost_ic =         l.out_scenario[["l.out_strategy"]][["int"]][["m.out"]][,c("cost_ic")], 
    cost    =         l.out_scenario[["l.out_strategy"]][["int"]][["m.out"]][,c("cost")], 
    nhb     =         l.out_scenario[["l.out_strategy"]][["int"]][["m.out"]][,c("nhb")]
  )
  a.trace <- array(data=c(m.trace_soc,m.trace_int), dim=c(nrow(m.trace_soc),ncol(m.trace_soc),2), dimnames=list(NULL,colnames(m.trace_soc),c("soc","int")))
  
  # initialize matrix for summary outcomes
  m.out_summary <- matrix(
    data = NA, 
    nrow = ncol(m.trace_soc), 
    ncol = 11, 
    dimnames = list(
      colnames(m.trace_soc),
      c("soc","int","dif","dif_relative","soc_within","int_within","soc_extrapolate","int_extrapolate","dif_within","dif_extrapolate","dif_p_extrapolate")
    )
  )
  m.out_summary
  
  # fill matrix
  for(i in colnames(m.trace_soc)) {
    m.out_summary[i,"soc"] <- sum(a.trace[,i,"soc"])
    m.out_summary[i,"int"] <- sum(a.trace[,i,"int"])
    m.out_summary[i,"soc_within"] <- sum(a.trace[1:2,i,"soc"])
    m.out_summary[i,"int_within"] <- sum(a.trace[1:2,i,"int"])
    m.out_summary[i,"soc_extrapolate"] <- sum(a.trace[3:29,i,"soc"])
    m.out_summary[i,"int_extrapolate"] <- sum(a.trace[3:29,i,"int"])
  }
  
  # add absolute, relative and proportional difference
  round(m.out_summary,1)
  m.out_summary[,"dif"] <- m.out_summary[,"int"] - m.out_summary[,"soc"]
  m.out_summary[,"dif_relative"] <- m.out_summary[,"dif"] / m.out_summary[,"soc"]
  m.out_summary[,"dif_within"] <- m.out_summary[,"int_within"] - m.out_summary[,"soc_within"]
  m.out_summary[,"dif_extrapolate"] <- m.out_summary[,"int_extrapolate"] - m.out_summary[,"soc_extrapolate"]
  m.out_summary[,"dif_p_extrapolate"] <- m.out_summary[,"dif_extrapolate"] / m.out_summary[,"dif"]
  round(m.out_summary, digits=2)
    # proportional difference is invalid for outcomes with negative values
  
  # short outcomes
  m.out_short <- matrix(data=NA, nrow=11, ncol=3, dimnames=list(c("mci","mil","mod","sev","alv","home","instit","ontx","cost","qaly","nhb"),c("soc","int","dif")))
  for(i in rownames(m.out_short)) {
    m.out_short[i,"soc"] <- sum(a.trace[1:n.cycles,i,"soc"])
    m.out_short[i,"int"] <- sum(a.trace[1:n.cycles,i,"int"])
  }
  m.out_short[,"dif"] <- m.out_short[,"int"] - m.out_short[,"soc"]
  
  # return
  return(list(
    a.trace = a.trace, 
    m.out_summary = m.out_summary, 
    m.out_short = m.out_short
  ))
  
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
    sex = "female", 
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
    p.discontinuation1 = 0, 
    p.discontinuation2 = 0.1, 
    discontinuation2_start = 2, 
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
    c.mci_sc = 222 * 12, 
    c.mil_sc = 410 * 12, 
    c.mod_sc = 653 * 12, 
    c.sev_sc = 1095 * 12, 
    c.mci_ic = 0, 
    c.mil_ic = 0, 
    c.mod_ic = 0, 
    c.sev_ic = 0, 
    c.mci_i_hc = 1254 * 12, 
    c.mil_i_hc = 1471 * 12, 
    c.mod_i_hc = 1958 * 12, 
    c.sev_i_hc = 2250 * 12, 
    c.mci_i_sc = 8762 * 12, 
    c.mil_i_sc = 8762 * 12, 
    c.mod_i_sc = 8762 * 12, 
    c.sev_i_sc = 8762 * 12, 
    c.mci_i_ic = 0, 
    c.mil_i_ic = 0, 
    c.mod_i_ic = 0, 
    c.sev_i_ic = 0, 
    c.Tx = 10000, 
    c.Tx_diagnostics1 = 2000, 
    discount_QALY = 0.035, 
    discount_COST = 0.035, 
    wtp = 40000, 
    half_cycle_correction = FALSE
  )
  
  # run model
  out_val_int <- f.run_scenario(l.inputs = l.inputs_val_int, detailed = TRUE)
  out_val_int
  
  ## incremental summary outcomes
  df.he_incr_val_int <- rbind(
    out_val_int[["df.out_sum"]][,-1], 
    incremental = out_val_int[["df.out_sum"]]["int",-1] - out_val_int[["df.out_sum"]]["soc",-1]
  )
  table.he_incr_val_int <- format(df.he_incr_val_int, digits=2, scientific=FALSE, big.mark=",")
  print(table.he_incr_val_int)
  
  # comment: fully replicated outcomes compared to previous version (main branch at this moment)

}



######################################## 4.1. EXTREME SCENARIOS ########################################
# to-do

######################################## 5. ANALYSIS ########################################

######################################## 5.1. PROBABILISTIC ########################################
# to be developed for this version, see release '2.1.0 ISPOR Europe 2023 abstract' for an older version


######################################## 5.2. DETERMINISTIC ########################################

######################################## 5.2.1. BASE CASE ########################################

# run the model
out_base <- f.run_scenario(l.inputs = l.inputs, detailed = TRUE)

# all outcomes
  # str(out_base) # show structure of output
  # out_base # all output
  # out_base[["l.out_strategy"]] # strategy details
  # out_base[["l.out_strategy"]][["soc"]] # strategy 'soc' details
  # out_base[["l.out_strategy"]][["int"]] # strategy 'int' details
  # out_base[["l.out_strategy"]][["soc"]][["a.TP"]] # strategy 'soc' transition probability matrix
  # out_base[["l.out_strategy"]][["soc"]][["m.trace"]] # strategy 'soc' state trace
  # out_base[["l.out_strategy"]][["soc"]][["m.out"]] # strategy 'soc' outcomes per cycle
  # out_base[["l.out_strategy"]][["int"]][["a.TP"]] # idem for 'int'
  # out_base[["l.out_strategy"]][["int"]][["m.trace"]]
  # out_base[["l.out_strategy"]][["int"]][["m.out"]]
  # out_base[["df.out_sum"]] # scenario results


round(out_base[["l.out_strategy"]][["soc"]][["m.trace"]],3)
round(colSums(out_base[["l.out_strategy"]][["soc"]][["m.trace"]]),2)
print(round(colSums(out_base[["l.out_strategy"]][["int"]][["m.trace"]]),2))

# summary outcomes
out_base_prepared <- f.prepare_outcomes(l.out_scenario = out_base, n.cycles = 29)
  print(round(out_base_prepared[["m.out_short"]],5))

# prepare: icer
icer <- calculate_icers(
  cost = out_base[["df.out_sum"]][,"COST"],
  effect = out_base[["df.out_sum"]][,"QALY"],
  strategies = out_base[["df.out_sum"]][,"strategy"]
)
print(icer)

# optional: copy outcomes to clipboard
  # temp.trace_soc <- out_base[["l.out_strategy"]][["soc"]][["m.trace"]]
  # write.table(x = temp.trace_soc, file = "clipboard", sep = "\t", row.names = FALSE)
  # temp.trace_int <- out_base[["l.out_strategy"]][["int"]][["m.trace"]]
  # write.table(x = temp.trace_int, file = "clipboard", sep = "\t", row.names = FALSE)


# print tables/plots

if(F) {
  
  ## table: summary outcomes
  table.out_summary <- format(out_base_prepared[["m.out_short"]], digits=2, scientific=FALSE, big.mark=",")
  print(table.out_summary)
  
  ## figure: state trace
  v.age_range <- c(l.inputs[["age_start"]]:(l.inputs[["age_start"]]+l.inputs[["n.cycle"]]-1)) # store age range
  xx <- c(v.age_range, rev(v.age_range)) # prepare polygon x-values
  a.trace <- out_base_prepared[["a.trace"]]
  yy_mci <- c(a.trace[,"mci","soc"], rev(a.trace[,"mci","int"])) # polygon y-values
  yy_mil <- c(a.trace[,"mil","soc"], rev(a.trace[,"mil","int"])) # idem
  yy_mod <- c(a.trace[,"mod","soc"], rev(a.trace[,"mod","int"])) # idem
  yy_sev <- c(a.trace[,"sev","soc"], rev(a.trace[,"sev","int"])) # idem
  yy_dth <- c(a.trace[,"dth","soc"], rev(a.trace[,"dth","int"])) # idem
  par(mar=c(5, 4, 4, 2)+0.1, xpd=FALSE)
  matplot(
    x = v.age_range, 
    y = a.trace[,,"soc"], 
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
    y = a.trace[,c("mci","mil","mod","sev","dth"),"soc"], 
    type = "l",
    lty = 1,
    col = c("green","yellow","orange","red","black")
  )
  matlines(
    x = v.age_range, 
    y = a.trace[,c("mci","mil","mod","sev","dth"),"int"], 
    type = "l",
    lty = 2,
    col = c("green","yellow","orange","red","black")
  )
  legend(x="topright", legend=c("mci","mil","mod","sev","dth"), col=c("green","yellow","orange","red","black"), lty=1)
  legend(x="right", legend=c("soc","int"), col="black", lty=c(1,2))
  
  ## figure: mean time in state
  m.time_state2 <- m.out_summary[c("mci","mil","mod","sev"),c("soc","int")]
  par(mar=c(8, 4, 4, 2), xpd=TRUE)
  barplot(
    height = m.time_state2, 
    horiz = TRUE, 
    xlab = "time (years)", 
    ylab = "strategy", 
    col=c("green","yellow","orange","red"), 
    space = 0.2, 
    main = "mean time in state"
  )
  legend(x="bottom", legend=c("mci","mil","mod","sev"), inset=c(0,-0.5), horiz=TRUE, fill=c("green","yellow","orange","red"))
  text(x=c(0,cumsum(m.time_state2[1:3,"soc"])), y=1, labels=round(m.time_state2[,"soc"],1), pos=4)
  text(x=c(0,cumsum(m.time_state2[1:3,"int"])), y=2, labels=round(m.time_state2[,"int"],1), pos=4)
  
  ## table: proportion in state
  print(round(a.trace[1:10,c("mci","mil","mod","sev"),"soc"], 2)) # state trace standard of care strategy
  print(round(a.trace[1:10,c("mci","mil","mod","sev"),"int"], 2)) # state trace intervention strategy
  
  ## figure: icer
  par(mar=c(5, 4, 4, 2)+0.1, xpd=FALSE)
  #print(plot(icer, label="all"))
  
  ## table: icer
  print(as.data.frame(t(icer)))
  
  ## figure: annual cost difference by sector over time
  m.cost_incr.pos <- m.cost_incr.neg <- m.trace_int[,c("cost_dx","cost_tx","cost_hc","cost_sc","cost_ic")] - m.trace_soc[,c("cost_dx","cost_tx","cost_hc","cost_sc","cost_ic")] # split positive and negative
  m.cost_incr.pos[m.cost_incr.pos<0] <- 0
  m.cost_incr.neg[m.cost_incr.neg>=0] <- 0
  par(mar=c(5, 4, 4, 2)+0.1, xpd=FALSE)
  barplot(
    height = t(m.cost_incr.pos),
    beside = F,
    xlab = "time (years)",
    ylab = "annual incremental costs",
    ylim = c(-5000, 6000),
    col = rainbow(5), 
    names.arg = c(NA,NA,NA,NA,5,NA,NA,NA,NA,10,NA,NA,NA,NA,15,NA,NA,NA,NA,20,NA,NA,NA,NA,25,NA,NA,NA,NA),
    main = "costs by sector over time"
  )
  barplot(
    height = t(m.cost_incr.neg),
    beside = F,
    col = rainbow(5),
    add = T
  )
  legend(x = "topright", legend = c("health","social","informal","treatment","diagnostic"), fill = rainbow(5))

}



######################################## 5.2.2. CYCLE TIME ########################################

if(F) {

  # cycle time adjustment following guidance by Gidwani et al. [2020: https://doi.org/10.1007/s40273-020-00937-z]
  
  # function to convert probability (p) to a different cycle time (t=time proportional to current 1-year cycle length)
  f.p_time <- function(p, t) 1-(1-p)^(t)
  
  # proportion of the 1-year cycle length
  t <- 6/12
  
  # simultaneously convert matrix of dementia transition probabilities to new cycle length
  ## temporary matrix of dementia transition probabilities
  t.TP <- matrix(data=NA, nrow=3, ncol=3, dimnames=list(c("mil","mod","sev"),c("mil","mod","sev")))
  t.TP["mil","mod"] <- l.inputs[["p.mil_mod"]]
  t.TP["mil","sev"] <- l.inputs[["p.mil_sev"]]
  t.TP["mod","mil"] <- l.inputs[["p.mod_mil"]]
  t.TP["mod","sev"] <- l.inputs[["p.mod_sev"]]
  t.TP["sev","mil"] <- l.inputs[["p.sev_mil"]]
  t.TP["sev","mod"] <- l.inputs[["p.sev_mod"]]
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
  t.TP2[t.TP2 < 0] <- 0 #convert small negative values to 0
  t.TP2 <- round(t.TP2, 3) # round
  dimnames(t.TP2) <- dimnames(t.TP) # add names
  rowSums(t.TP2) # check if rowsums add up to 1
  diag(t.TP2) <- NA # remove matrix diagonal (probability of remaining in the same state)
  t.TP2[1,1] <- 1-sum(t.TP2[1,-1]) # calculate probability of remaining in the same state as '1-(other probabilities)'
  t.TP2[2,2] <- 1-sum(t.TP2[2,-2])
  t.TP2[3,3] <- 1-sum(t.TP2[3,-3])
  t.TP2
  
  # check: compare p.mci_mil handled separately versus handled simultaneously with dementia transition probabilities
  ## handled separately
  t.TP2_mci <- rbind(0,cbind(0,t.TP2))
  rownames(t.TP2_mci)[1] <- "mci"
  colnames(t.TP2_mci)[1] <- "mci"
  t.TP2_mci["mci","mil"] <- f.p_time(p = l.inputs[["p.mci_mil"]], t = t)
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
  t.TP_mci["mci","mil"] <- l.inputs[["p.mci_mil"]]
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
  
  # adjust mortality table
  m.mortality_rate_US_cycle4 <- m.mortality_rate_US * t # mortality table is rate, which can be multiplied rather than 
  t.rows <- rep(x = 1:nrow(m.mortality_rate_US_cycle4), each = 1/t) # create vector for row numbers and repeat each row multiple times
  m.mortality_rate_US_cycle4 <- m.mortality_rate_US_cycle4[t.rows,] # select each row multiple times
  
  # adjust inputs to cycle time (list is duplicated to force check all parameters for adjustment)
  l.inputs_c4 <- list(
    v.names_state = c("mcion","mciof","milon","milof","mod","sev","mci_i","mil_i","mod_i","sev_i","dth"), 
    v.names_strat = c("soc","int"), 
    age_start = 70/t, 
    n.cycle = 29/t, 
    sex = "female", 
    p.mci_mil = f.p_time(p = l.inputs[["p.mci_mil"]], t = t), 
    p.mci_mod = 0, # if this probability is different from 0, we think it is better to include it in the simultaneous conversion of transition probabilities to cycle time (as done for dementia transitions above)
    p.mci_sev = 0, # idem
    p.mil_mci = 0, # idem
    p.mil_mod = t.TP2["mil","mod"], 
    p.mil_sev = t.TP2["mil","sev"], 
    p.mod_mil = t.TP2["mod","mil"], 
    p.mod_sev = t.TP2["mod","sev"], 
    p.sev_mil = t.TP2["sev","mil"], 
    p.sev_mod = t.TP2["sev","mod"], 
    p.mci_i = f.p_time(p = l.inputs[["p.mci_i"]], t = t), 
    p.mil_i = f.p_time(p = l.inputs[["p.mil_i"]], t = t), 
    p.mod_i = f.p_time(p = l.inputs[["p.mod_i"]], t = t), 
    p.sev_i = f.p_time(p = l.inputs[["p.sev_i"]], t = t), 
    p.discontinuation1 = f.p_time(p = l.inputs[["p.discontinuation1"]], t = t), 
    p.discontinuation2 = f.p_time(p = l.inputs[["p.discontinuation2"]], t = t), 
    discontinuation2_start = l.inputs[["discontinuation2_start"]], 
    m.r.mortality = m.mortality_rate_US_cycle4, 
    hr.mort_mci = 1, 
    hr.mort_verymilddem = 1.82, 
    hr.mort_mil = 1.318, 
    hr.mort_mod = 2.419, 
    hr.mort_sev = 4.267, 
    rr.tx_mci_mil = 0.75, 
    rr.tx_mci_mod = 1, 
    rr.tx_mci_sev = 1, 
    rr.tx_mil_mod = 0.75, 
    rr.tx_mci_mil_dis = 0.75, 
    rr.tx_mci_mod_dis = 1, 
    rr.tx_mci_sev_dis = 1, 
    rr.tx_mil_mod_dis = 0.75, 
    tx_waning = f.p_time(p = l.inputs[["tx_waning"]], t = t), 
    tx_waning_dis = f.p_time(p = l.inputs[["tx_waning_dis"]], t = t), 
    tx_duration = l.inputs[["tx_duration"]] * (1/t), 
    p.starting_state_mci = 1, 
    u.mci_pt = l.inputs[["u.mci_pt"]] * t, # adjusted as summing of QALYs is not adjusted for cycle length in the model code
    u.mil_pt = l.inputs[["u.mil_pt"]] * t, # idem
    u.mod_pt = l.inputs[["u.mod_pt"]] * t, # idem
    u.sev_pt = l.inputs[["u.sev_pt"]] * t, # idem
    u.mci_ic = l.inputs[["u.mci_ic"]] * t, # idem
    u.mil_ic = l.inputs[["u.mil_ic"]] * t, # idem
    u.mod_ic = l.inputs[["u.mod_ic"]] * t, # idem
    u.sev_ic = l.inputs[["u.sev_ic"]] * t, # idem
    c.mci_hc = l.inputs[["c.mci_hc"]] * t, 
    c.mil_hc = l.inputs[["c.mil_hc"]] * t, 
    c.mod_hc = l.inputs[["c.mod_hc"]] * t, 
    c.sev_hc = l.inputs[["c.sev_hc"]] * t, 
    c.mci_sc = l.inputs[["c.mci_sc"]] * t, 
    c.mil_sc = l.inputs[["c.mil_sc"]] * t, 
    c.mod_sc = l.inputs[["c.mod_sc"]] * t, 
    c.sev_sc = l.inputs[["c.sev_sc"]] * t, 
    c.mci_ic = l.inputs[["c.mci_ic"]] * t, 
    c.mil_ic = l.inputs[["c.mil_ic"]] * t, 
    c.mod_ic = l.inputs[["c.mod_ic"]] * t, 
    c.sev_ic = l.inputs[["c.sev_ic"]] * t, 
    c.mci_i_hc = l.inputs[["c.mci_hc"]] * t, 
    c.mil_i_hc = l.inputs[["c.mil_hc"]] * t, 
    c.mod_i_hc = l.inputs[["c.mod_hc"]] * t, 
    c.sev_i_hc = l.inputs[["c.sev_hc"]] * t, 
    c.mci_i_sc = l.inputs[["c.mci_sc"]] * t, 
    c.mil_i_sc = l.inputs[["c.mil_sc"]] * t, 
    c.mod_i_sc = l.inputs[["c.mod_sc"]] * t, 
    c.sev_i_sc = l.inputs[["c.sev_sc"]] * t, 
    c.mci_i_ic = l.inputs[["c.mci_ic"]] * t, 
    c.mil_i_ic = l.inputs[["c.mil_ic"]] * t, 
    c.mod_i_ic = l.inputs[["c.mod_ic"]] * t, 
    c.sev_i_ic = l.inputs[["c.sev_ic"]] * t, 
    c.Tx = l.inputs[["c.Tx"]] * t, 
    c.Tx_diagnostics1 = 2000, 
    discount_QALY = f.p_time(p = l.inputs[["discount_QALY"]], t = t), 
    discount_COST = f.p_time(p = l.inputs[["discount_COST"]], t = t), 
    wtp = 40000 
  )
  
  # run the model
  out_c4 <- f.run_scenario(l.inputs = l.inputs_c4, detailed = TRUE)
  
  # compare state trace
  round(out_base[["l.out_strategy"]][["soc"]][["m.trace"]][c(1:5,20),],2)
  round(out_c4  [["l.out_strategy"]][["soc"]][["m.trace"]][c(1:5,20)*(1/t)-1,],2)
  # compare health-economic outcomes
  out_base[["df.out_sum"]]
  out_c4[["df.out_sum"]] # !!!TO-DO: life years must be corrected for different cycle time
  # !!!TO-DO: these results need to be further validated

}


######################################## 5.2.3. DETERMINISTIC SENSITIVITY ANALYSIS ########################################

if(F) {
  
  # list parameters for DSA
  dsa_pars <- c(
    "age_start",
    "hr.mort_verymilddem",
    "hr.mort_mil",
    "hr.mort_mod",
    "hr.mort_sev",
    "p.mci_mil",
    "rr.tx_mci_mil",
    "rr.tx_mil_mod",
    "tx_duration",
    "c.Tx"
  )
  
  # list minimum values
  dsa_min <- c(
    60,
    1.55,
    1.153,
    2.122,
    3.610,
    l.inputs[["p.mci_mil"]]/2,
    l.inputs[["rr.tx_mci_mil"]]^2,
    l.inputs[["rr.tx_mil_mod"]]^2,
    l.inputs[["tx_duration"]]/2,
    l.inputs[["c.Tx"]]/5
  )
  
  # list maximum values
  dsa_max <- c(
    80,
    2.14,
    1.507,
    2.757,
    5.043,
    l.inputs[["p.mci_mil"]]*2,
    l.inputs[["rr.tx_mci_mil"]]^0.5,
    l.inputs[["rr.tx_mil_mod"]]^0.5,
    l.inputs[["tx_duration"]]*2,
    l.inputs[["c.Tx"]]*5
  )
  
  # define sensitivity analysis range
  df.owsa_params_range <- data.frame(
    pars = dsa_pars,
    min = dsa_min,
    max = dsa_max
  ); df.owsa_params_range
  
  # run DSA
  out_owsa_det <- run_owsa_det(
    params_range = df.owsa_params_range,
    params_basecase = l.inputs,
    nsamp = 25, # n equally-spaced samples
    FUN = f.run_scenario, # make sure the f.run_scenario function only returns a single dataframe with aggregated outcomes (not detailed outcomes)
    outcomes = NULL, # NULL for all outcomes (requires to select specific outcome for plotting the tornado graph later)
    strategies = NULL, # NULL for all strategies
    progress = TRUE
  )
  str(out_owsa_det) # show structure of result
  out_owsa_det # show results
  
  # !! IMPORTANT: the function 'owsa_tornado' omits parameters with negative values due to the 'min_rel_diff' option. Therefore, you need to manually adjust the function. The following is suggested: 
  # step 1: run the following code: trace("owsa_tornado", edit=TRUE)
  # step 2: manually adjust the code line 20 by adding 'abs()' function: 'abs(.data$outcome_val.low)' instead of '.data$outcome_val.low'
  
  # plot tornado graph of DSA outcome 'NMB' (fyi: not incremental)
  select_owsa <- subset(x=out_owsa_det$owsa_NHB, subset=strategy=="st_int") # select outcome and from which strategy
  owsa_tornado(select_owsa)
  
  # plot tornado graph of incremental outcome (only use when there are 2 strategies, in case of more than 2 strategies see 'dampack' vignette for notes on tornado limitation: https://cran.r-project.org/web/packages/dampack/)
  out_owsa_det_NHB_soc <- subset(x=out_owsa_det$owsa_NHB, subset=strategy=="st_soc") # select outcome and from which strategy
  out_owsa_det_NHB_int <- subset(x=out_owsa_det$owsa_NHB, subset=strategy=="st_int") # select outcome and from alternative strategy
  if(!all(out_owsa_det_NHB_soc[,c("parameter","param_val")]==out_owsa_det_NHB_int[,c("parameter","param_val")])) stop("error") # check to determine if data is correctly compared
  out_owsa_det_iNHB <- out_owsa_det_NHB_int
  out_owsa_det_iNHB[,"outcome_val"] <- NA
  out_owsa_det_iNHB[,"outcome_val"] <- out_owsa_det_NHB_int[,"outcome_val"] - out_owsa_det_NHB_soc[,"outcome_val"]
  owsa_tornado(out_owsa_det_iNHB, n_y_ticks=4, return = "data")
  owsa_tornado(out_owsa_det_iNHB, n_y_ticks=4, return = "plot")
  
  # plot optimal strategy (outcome over each parameter range for each strategy)
  plot(out_owsa_det$owsa_NHB, n_x_ticks = 3)
  # Visualize optimal strategy (max NMB) over each parameter range
  owsa_opt_strat(out_owsa_det$owsa_NHB, return="data")
  #windows()
  owsa_opt_strat(out_owsa_det$owsa_NHB, return="plot")
  
}



# [to be implemented in updated version]


# natural progression MCI to dementia
## source: [Vos, 2015: https://doi.org/10.1093/brain/awv029]

## Amyloid positive & neuronal loss positive
### Operationalized by diagnostic criteria NIA-AA categories: 'NIA-AA high AD' (Amyloid+, Injury+)
### corresponding 3-year cumulative incidence probability: 'high AD' = 59% (AD dementia; table 3) and 4% (non-AD dementia; mentioned in text)
### This results into a weighted 3-year cumulative incidence of:
temp.est1 <- 1-exp(- (-log(1-0.59) + -log(1-0.04)) )
## and corresponding 1-year probability of 
temp.est1 <- 1-exp(- -log(1-temp.est1)/3)
temp.est1

## Amyloid positive & neuronal loss undetermined
### Operationalized by diagnostic criteria NIA-AA categories: 'NIA-AA high AD' (Amyloid+, Injury+) and 'conflicting IAP' (Amyloid+, Injury-)
### corresponding 3-year cumulative incidence probability: 'high AD' = 59% (AD dementia; table 3) and 4% (non-AD dementia; mentioned in text), and 22% (AD dementia; table 3) and 4% (non-AD dementia; mentioned in text) with prevalence of 353 and 49 respectively (respectively)
### This results into a weighted 3-year cumulative incidence of (converting all 4 probabilities to rates before weighting and averaging):
temp.est2 <- 1-exp(- ( (-log(1-0.59) + -log(1-0.04))*353 + (-log(1-0.22) + -log(1-0.04))*49 ) / (353+49) )
## and corresponding 1-year probability of 
temp.est2 <- 1-exp(- -log(1-temp.est2)/3)
temp.est2




######################################## 5.2.4. HEADROOM ########################################
# to be developed for this version, see release '2.1.0 ISPOR Europe 2023 abstract' for an older version



######################################## 5.2.5. ICER REPLICATION ########################################

if(T) {
  
  # input parameters
  l.inputs_icer <- list(
    v.names_state = c("mcion","mciof","milon","milof","mod","sev","mci_i","mil_i","mod_i","sev_i","dth"), # disease states: mci = mild cognitive impairment; mil = mild dementia; mod = moderate dementia; sev = severe dementia; dth = dead; x_i = living in institutional setting (without '_i' = living in community)
    v.names_strat = c("soc","int"), # strategies: soc = standard of care strategy; int = intervention strategy
    age_start = 70, 
    n.cycle = 29, 
    sex = "male", 
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
    m.r.mortality = m.mortality_rate_US, 
    hr.mort_mci = 1.82, 
    hr.mort_verymilddem = 1, 
    hr.mort_mil = 2.92, 
    hr.mort_mod = 3.85, 
    hr.mort_sev = 9.52, 
    rr.tx_mci_mil = 0.69, 
    rr.tx_mci_mod = 1, 
    rr.tx_mci_sev = 1, 
    rr.tx_mil_mod = 0.69, 
    rr.tx_mci_mil_dis = 0.69, 
    rr.tx_mci_mod_dis = 1, 
    rr.tx_mci_sev_dis = 1, 
    rr.tx_mil_mod_dis = 0.69, 
    p.discontinuation1 = 0.069, 
    p.discontinuation2 = 0, 
    discontinuation2_start = 2, 
    tx_waning = 0, 
    tx_waning_dis = 0, 
    tx_duration = 29, 
    p.starting_state_mci = 0.55, 
    u.mci_pt = 0.73, 
    u.mil_pt = 0.73 + 0.17 - 0.22, 
    u.mod_pt = 0.73 + 0.17 - 0.36, 
    u.sev_pt = 0.73 + 0.17 - 0.53, 
    u.mci_ic = -0.03, 
    u.mil_ic = -0.05, 
    u.mod_ic = -0.08, 
    u.sev_ic = -0.10, 
    c.mci_hc = 1254 * 12, 
    c.mil_hc = 1471 * 12, 
    c.mod_hc = 1958 * 12, 
    c.sev_hc = 2250 * 12, 
    c.mci_sc = 222 * 12, 
    c.mil_sc = 410 * 12, 
    c.mod_sc = 653 * 12, 
    c.sev_sc = 1095 * 12, 
    c.mci_ic =  988 * 12, 
    c.mil_ic = 2184 * 12, 
    c.mod_ic = 3227 * 12, 
    c.sev_ic = 5402 * 12, 
    c.mci_i_hc = 1254 * 12, 
    c.mil_i_hc = 1471 * 12, 
    c.mod_i_hc = 1958 * 12, 
    c.sev_i_hc = 2250 * 12, 
    c.mci_i_sc = 8762 * 12, 
    c.mil_i_sc = 8762 * 12, 
    c.mod_i_sc = 8762 * 12, 
    c.sev_i_sc = 8762 * 12, 
    c.mci_i_ic =  435 * 12, 
    c.mil_i_ic =  961 * 12, 
    c.mod_i_ic = 1420 * 12, 
    c.sev_i_ic = 2377 * 12, 
    c.Tx = 26500 + (52/2)*78.35, # drug annual wholesale acquisition cost + treatment administration frequency * administration cost
    c.Tx_diagnostics1 = 261.10*3 + 261.10*3*0.215, # mri cost * 3-month monitoring in year 1 + mri cost * 3 times * proportion aria
    discount_EFFECT = 0.03, 
    discount_QALY = 0.03, 
    discount_COST = 0.03, 
    wtp = 100000, 
    half_cycle_correction = FALSE
  )
  
  # run the model
  out_base <- f.run_scenario(l.inputs = l.inputs_icer, detailed = TRUE)
  
  # discount life years
  ly_soc <- out_base[["l.out_strategy"]][["soc"]][["m.out"]][,"ly"]
  ly_int <- out_base[["l.out_strategy"]][["int"]][["m.out"]][,"ly"]
  v.discount <- 1 / (( 1 + (0.03)) ^ (0:28))
  sum(v.discount * ly_soc)
  sum(v.discount * ly_int)
  
  sum(v.discount * ly_int) - sum(v.discount * ly_soc)
  
 
}


######################################## 5.2.6. PREPARE SENSITIVITY ANALYSIS ########################################

# s_cdrhr
l.inputs_s_cdrhr <- l.inputs_icer
l.inputs_s_cdrhr[["rr.tx_mci_mil"]] = 0.66
l.inputs_s_cdrhr[["rr.tx_mil_mod"]] = 0.66
l.inputs_s_cdrhr[["rr.tx_mci_mil_dis"]] = 0.66
l.inputs_s_cdrhr[["rr.tx_mil_mod_dis"]] = 0.66
l.inputs_s_cdrhr[["p.discontinuation1"]] = 0.069
l.inputs_s_cdrhr[["p.discontinuation2"]] = 0
l.inputs_s_cdrhr[["discontinuation2_start"]] = 2
l.inputs_s_cdrhr[["tx_waning"]] = 0
l.inputs_s_cdrhr[["tx_waning_dis"]] = 0
l.inputs_s_cdrhr[["tx_duration"]] = 29
l.inputs_s_cdrhr[["p.starting_state_mci"]] = 0.55

# run model
out_s_cdrhr <- f.run_scenario(l.inputs = l.inputs_s_cdrhr, detailed = TRUE)

# run results
out_s_cdrhr_prepared <- f.prepare_outcomes(l.out_scenario = out_s_cdrhr, n.cycles = 29)

# s_cdr_%
# s_cdr_timeshift
# s_mmse_hr
# s_mmse_%
# s_mmse_time

# sA (= cdrhr)
l.inputs_sA <- l.inputs_s_cdrhr
l.inputs_sA[["rr.tx_mci_mil_dis"]] = 0.66
l.inputs_sA[["rr.tx_mil_mod_dis"]] = 0.66
l.inputs_sA[["p.discontinuation2"]] = 0
l.inputs_sA[["discontinuation2_start"]] = 2
l.inputs_sA[["tx_waning"]] = 0
l.inputs_sA[["tx_waning_dis"]] = 0
out_sA <- f.run_scenario(l.inputs = l.inputs_sA, detailed = TRUE)
out_sA_prepared <- f.prepare_outcomes(l.out_scenario = out_sA, n.cycles = 29)

# sB
round((1-0.25)^c(0:7),2) # reflects decrease to relatively small treatment effect over 7 years
l.inputs_sB <- l.inputs_s_cdrhr
l.inputs_sB[["rr.tx_mci_mil_dis"]] = 0.66
l.inputs_sB[["rr.tx_mil_mod_dis"]] = 0.66
l.inputs_sB[["p.discontinuation2"]] = 0
l.inputs_sB[["discontinuation2_start"]] = 2
l.inputs_sB[["tx_waning"]] = 0.25
l.inputs_sB[["tx_waning_dis"]] = 0.25
l.inputs_sB[["tx_duration"]] = 29
out_sB <- f.run_scenario(l.inputs = l.inputs_sB, detailed = TRUE)
out_sB_prepared <- f.prepare_outcomes(l.out_scenario = out_sB, n.cycles = 29)

# sC
l.inputs_sC <- l.inputs_s_cdrhr
l.inputs_sC[["rr.tx_mci_mil_dis"]] = 1
l.inputs_sC[["rr.tx_mil_mod_dis"]] = 1
l.inputs_sC[["p.discontinuation2"]] = 1
l.inputs_sC[["discontinuation2_start"]] = 2
l.inputs_sC[["tx_waning"]] = 0
l.inputs_sC[["tx_waning_dis"]] = 0
l.inputs_sC[["tx_duration"]] = 29
out_sC <- f.run_scenario(l.inputs = l.inputs_sC, detailed = TRUE)
out_sC_prepared <- f.prepare_outcomes(l.out_scenario = out_sC, n.cycles = 29)

# sD
l.inputs_sD <- l.inputs_s_cdrhr
l.inputs_sD[["rr.tx_mci_mil_dis"]] = 0.66
l.inputs_sD[["rr.tx_mil_mod_dis"]] = 0.66
l.inputs_sD[["p.discontinuation2"]] = 1
l.inputs_sD[["discontinuation2_start"]] = 2
l.inputs_sD[["tx_waning"]] = 0
l.inputs_sD[["tx_waning_dis"]] = 0
l.inputs_sD[["tx_duration"]] = 29
out_sD <- f.run_scenario(l.inputs = l.inputs_sD, detailed = TRUE)
out_sD_prepared <- f.prepare_outcomes(l.out_scenario = out_sD, n.cycles = 29)

# sE
l.inputs_sE <- l.inputs_s_cdrhr
l.inputs_sE[["rr.tx_mci_mil_dis"]] = 0.66
l.inputs_sE[["rr.tx_mil_mod_dis"]] = 0.66
l.inputs_sE[["p.discontinuation2"]] = 1
l.inputs_sE[["discontinuation2_start"]] = 2
l.inputs_sE[["tx_waning"]] = 0.25
l.inputs_sE[["tx_waning_dis"]] = 0.25
l.inputs_sE[["tx_duration"]] = 29
out_sE <- f.run_scenario(l.inputs = l.inputs_sE, detailed = TRUE)
out_sE_prepared <- f.prepare_outcomes(l.out_scenario = out_sE, n.cycles = 29)

# initialize: summary incremental outcomes
m.sABCDE <- matrix(data = NA, nrow = 5, ncol = 6, dimnames = list(c("A","B","C","D","E"),c("mcimil","alv","qaly_pt","cost_dxtx","cost_hcscic","nhb")))
# fill
m.sABCDE["A","mcimil"]      <- sum(out_sA_prepared[["m.out_summary"]][c("mci","mil"),"dif"])
m.sABCDE["A","alv"]         <-     out_sA_prepared[["m.out_summary"]]["alv","dif"]
m.sABCDE["A","qaly_pt"]     <-     out_sA_prepared[["m.out_summary"]]["qaly_pt","dif"]
m.sABCDE["A","cost_dxtx"]   <- sum(out_sA_prepared[["m.out_summary"]][c("cost_dx","cost_tx"),"dif"])
m.sABCDE["A","cost_hcscic"] <- sum(out_sA_prepared[["m.out_summary"]][c("cost_hc","cost_sc","cost_ic"),"dif"])
m.sABCDE["A","nhb"]         <-     out_sA_prepared[["m.out_summary"]][c("nhb"),"dif"]

m.sABCDE["B","mcimil"]      <- sum(out_sB_prepared[["m.out_summary"]][c("mci","mil"),"dif"])
m.sABCDE["B","alv"]         <-     out_sB_prepared[["m.out_summary"]]["alv","dif"]
m.sABCDE["B","qaly_pt"]     <-     out_sB_prepared[["m.out_summary"]]["qaly_pt","dif"]
m.sABCDE["B","cost_dxtx"]   <- sum(out_sB_prepared[["m.out_summary"]][c("cost_dx","cost_tx"),"dif"])
m.sABCDE["B","cost_hcscic"] <- sum(out_sB_prepared[["m.out_summary"]][c("cost_hc","cost_sc","cost_ic"),"dif"])
m.sABCDE["B","nhb"]         <-     out_sB_prepared[["m.out_summary"]][c("nhb"),"dif"]

m.sABCDE["C","mcimil"]      <- sum(out_sC_prepared[["m.out_summary"]][c("mci","mil"),"dif"])
m.sABCDE["C","alv"]         <-     out_sC_prepared[["m.out_summary"]]["alv","dif"]
m.sABCDE["C","qaly_pt"]     <-     out_sC_prepared[["m.out_summary"]]["qaly_pt","dif"]
m.sABCDE["C","cost_dxtx"]   <- sum(out_sC_prepared[["m.out_summary"]][c("cost_dx","cost_tx"),"dif"])
m.sABCDE["C","cost_hcscic"] <- sum(out_sC_prepared[["m.out_summary"]][c("cost_hc","cost_sc","cost_ic"),"dif"])
m.sABCDE["C","nhb"]         <-     out_sC_prepared[["m.out_summary"]][c("nhb"),"dif"]

m.sABCDE["D","mcimil"]      <- sum(out_sD_prepared[["m.out_summary"]][c("mci","mil"),"dif"])
m.sABCDE["D","alv"]         <-     out_sD_prepared[["m.out_summary"]]["alv","dif"]
m.sABCDE["D","qaly_pt"]     <-     out_sD_prepared[["m.out_summary"]]["qaly_pt","dif"]
m.sABCDE["D","cost_dxtx"]   <- sum(out_sD_prepared[["m.out_summary"]][c("cost_dx","cost_tx"),"dif"])
m.sABCDE["D","cost_hcscic"] <- sum(out_sD_prepared[["m.out_summary"]][c("cost_hc","cost_sc","cost_ic"),"dif"])
m.sABCDE["D","nhb"]         <-     out_sD_prepared[["m.out_summary"]][c("nhb"),"dif"]

m.sABCDE["E","mcimil"]      <- sum(out_sE_prepared[["m.out_summary"]][c("mci","mil"),"dif"])
m.sABCDE["E","alv"]         <-     out_sE_prepared[["m.out_summary"]]["alv","dif"]
m.sABCDE["E","qaly_pt"]     <-     out_sE_prepared[["m.out_summary"]]["qaly_pt","dif"]
m.sABCDE["E","cost_dxtx"]   <- sum(out_sE_prepared[["m.out_summary"]][c("cost_dx","cost_tx"),"dif"])
m.sABCDE["E","cost_hcscic"] <- sum(out_sE_prepared[["m.out_summary"]][c("cost_hc","cost_sc","cost_ic"),"dif"])
m.sABCDE["E","nhb"]         <-     out_sE_prepared[["m.out_summary"]][c("nhb"),"dif"]

print(round(m.sABCDE,2))

