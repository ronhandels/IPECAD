
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

# see readme.md for a short explanation of each input parameter

l.inputs <- list(
  v.names_state = c("mcion","mciof","milon","milof","mod","sev","mci_i","mil_i","mod_i","sev_i","dth"), # disease states: mci = mild cognitive impairment; mil = mild dementia; mod = moderate dementia; sev = severe dementia; dth = dead; x_i = living in institutional setting (without '_i' = living in community)
  v.names_strat = c("soc","int"), # strategies: soc = standard of care strategy; int = intervention strategy
  age_start = 70, 
  age_end = 99, 
  sex = "female", 
  p.mci_mil = 0.206, 
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
  p.discontinuation1 = 0.1, 
  p.discontinuation_x = 0.1, 
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
  tx_waning = 0.05, 
  tx_duration = 7, 
  p.starting_state_mci = 1, 
  u.mci = 0.73, 
  u.mil = 0.69, 
  u.mod = 0.53, 
  u.sev = 0.38, 
  c.mci = (1254 +  222) * 12, 
  c.mil = (1471 +  410) * 12, 
  c.mod = (1958 +  653) * 12, 
  c.sev = (2250 + 1095) * 12, 
  c.mci_i = (1254 + 8762) * 12, 
  c.mil_i = (1471 + 8762) * 12, 
  c.mod_i = (1958 + 8762) * 12, 
  c.sev_i = (2250 + 8762) * 12, 
  c.Tx = 10000, 
  c.Tx_diagnostics1 = 2000, 
  discount_QALY = 0.035, 
  discount_COST = 0.035, 
  wtp = 40000 
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
    
    # initialize output table (STEP G3: initialize objects to store strategy outcomes)
    m.out <- matrix(data = NA, nrow = n.cycle, ncol = 4, dimnames = list(NULL, c("qaly","cost","ly","nhb")))
    
    # initialize state trace
    m.trace <- matrix(data = NA, nrow = n.cycle, ncol = n.state, dimnames = list(NULL,v.names_state))
    
    # set starting state distribution (STEP G4: starting state)
    m.trace[1,] <- m.trace1
    
    # markov multiplication (STEP G5: markov multiplication by looping over cycles)
    for(t in 1:(n.cycle-1)) {
      m.trace[t+1,] <- m.trace[t,] %*% a.TP[,,t]
    }
    
    # primary economic outputs (STEP G6: multiply states with utility and cost estimates)
    m.out[,"ly"]   <- m.trace %*% c(1           , 1    , 1           , 1    , 1    , 1    , 1      , 1      , 1      , 1      , 0) # must match order of states
    m.out[,"qaly"] <- m.trace %*% c(u.mci       , u.mci, u.mil       , u.mil, u.mod, u.sev, u.mci  , u.mil  , u.mod  , u.sev  , 0) # must match order of states
    m.out[,"cost"] <- m.trace %*% c(c.mci + c.Tx, c.mci, c.mil + c.Tx, c.mil, c.mod, c.sev, c.mci_i, c.mil_i, c.mod_i, c.sev_i, 0) # must match order of states
    
    # half-cycle correction (STEP G7: apply half-cycle correction)
    for (i in 1:(n.cycle-1)) {
      m.out[i,"qaly"] <- (m.out[i,"qaly"] + m.out[i+1,"qaly"]) * 0.5
      m.out[i,"cost"] <- (m.out[i,"cost"] + m.out[i+1,"cost"]) * 0.5
      m.out[i,"ly"]   <- (m.out[i,"ly"]   + m.out[i+1,"ly"])   * 0.5
      m.out[i,"nhb"]  <- (m.out[i,"nhb"]  + m.out[i+1,"nhb"])  * 0.5
    }
    m.out[n.cycle,] <- 0
    
    # add additional diagnostics costs
    if(strat=="int") {
      m.out[,"cost"][1] <- m.out[,"cost"][1] + c.Tx_diagnostics1
    }
    
    # define vector for discounting QALYs and costs (STEP G8: apply discounting)
    v.discount_QALY <- 1 / (( 1 + (discount_QALY)) ^ (0 : (age_end-age_start-1)))
    v.discount_COST <- 1 / (( 1 + (discount_COST)) ^ (0 : (age_end-age_start-1)))
    
    # apply discounting
    m.out[,"qaly"] <- m.out[,"qaly"]*v.discount_QALY
    m.out[,"cost"] <- m.out[,"cost"]*v.discount_COST
    
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
    n.cycle <- (age_end-age_start) # number of cycles
    
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
      v.p.discontinuation[1] <- 0 # discontinuation at cycle 1 (starting states) fixed to 0
      v.p.discontinuation[2] <- p.discontinuation1 # discontinuation at cycle 2
      v.p.discontinuation[3:n.cycle] <- p.discontinuation_x # discontinuation at cycle 3 onwards
      v.p.discontinuation[tx_duration:n.cycle] <- 1 # maximum treatment duration implemented as discontinuation
      
      # death (subset mortality table to obtain age- and sex-specific mortality)
      v.r.dth <- m.r.mortality[age_start:(age_end-1), sex]
      
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
        
        # update transition probabilities treatment effect
        v.p.mcion_mil <- 1-exp(-(-log(1-p.mci_mil) * temp.rr.tx_mci_mil)) # convert probability to rate, then multiply with treatment relative risk, then convert to probability
        v.p.mcion_mod <- 1-exp(-(-log(1-p.mci_mod) * temp.rr.tx_mci_mod)) # idem
        v.p.mcion_sev <- 1-exp(-(-log(1-p.mci_sev) * temp.rr.tx_mci_sev)) # idem
        v.p.milon_mod <- 1-exp(-(-log(1-p.mil_mod) * temp.rr.tx_mil_mod)) # idem
        
        # update transition probabilities of remaining in the same state
        v.p.mcion_mci <- 1 - v.p.mcion_mil - v.p.mcion_mod - v.p.mcion_sev
        v.p.milon_mil <- 1 - v.p.milon_mci - v.p.milon_mod - v.p.mil_sev
        
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
    
    # return result
    if(!detailed) return(df.out_sum)
    if(detailed) {
      return(list(
        l.out_strategy = l.out_strategy,
        df.out_sum = df.out_sum
      ))
    }
  }
  )
}



######################################## 3. MODEL CALIBRATION ########################################
# n/a



######################################## 4. VALIDATION ########################################
# n/a



######################################## 5. ANALYSIS ########################################

######################################## 5.1. PROBABILISTIC ########################################
# to be developed for this version, see release '2.1.0 ISPOR Europe 2023 abstract' for an older version


######################################## 5.2. DETERMINISTIC ########################################


######################################## 5.2.1 BASE CASE ########################################

# run the model
out_base <- f.run_scenario(l.inputs = l.inputs, detailed = TRUE)

# possible output to choose from
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

# prepare data for tables/plots

## table: summary outcomes
temp.table1 <- rbind(
  out_base[["df.out_sum"]][,-1], 
  incremental=out_base[["df.out_sum"]]["int",-1] - out_base[["df.out_sum"]]["soc",-1] # calculate difference between 'soc' and 'int' strategies
)
df.table1 <- temp.table1
df.table1 <- format(temp.table1, digits=2, scientific=FALSE, big.mark=",") # format 

## plot: state trace
m.plot1_soc <- cbind(
  mci=out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"mcion"] + out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"mciof"] + out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"mci_i"], 
  mil=out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"milon"] + out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"milof"] + out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"mil_i"], 
  mod=out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"mod"] + out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"mod_i"], 
  sev=out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"sev"] + out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"sev_i"], 
  dth=out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"dth"]
)
m.plot1_int <- cbind(
  mci=out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"mcion"] + out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"mciof"] + out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"mci_i"], 
  mil=out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"milon"] + out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"milof"] + out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"mil_i"], 
  mod=out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"mod"] + out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"mod_i"], 
  sev=out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"sev"] + out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"sev_i"], 
  dth=out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"dth"]
)
a.plot1 <- array(data=c(m.plot1_soc,m.plot1_int), dim=c(nrow(m.plot1_soc),5,2), dimnames=list(NULL,colnames(m.plot1_soc),c("soc","int")))

## plot: mean time in state
m.plot2 <- cbind(soc=colSums(m.plot1_soc), int=colSums(m.plot1_int))
m.plot2 <- m.plot2[c("mci","mil","mod","sev"),]
# temp.dif <- t(t(m.plot2[,"int"] - m.plot2[,"soc"])) # difference in time in state (as column)
# m.plot2_neg <- m.plot2_pos <- temp.dif # split positive and negative differences
# m.plot2_neg[m.plot2_neg>=0] <- 0
# m.plot2_pos[m.plot2_pos<0 ] <- 0

## plot: icer
icer <- calculate_icers(
  cost = out_base[["df.out_sum"]][,"COST"],
  effect = out_base[["df.out_sum"]][,"QALY"],
  strategies = out_base[["df.out_sum"]][,"strategy"]
)


# print tables/plots

## table: summary outcomes
print(df.table1)

# plot: state trace
v.age_range <- c(l.inputs[["age_start"]]:(l.inputs[["age_end"]]-1)) # store age range
xx <- c(v.age_range, rev(v.age_range)) # prepare polygon x-values
yy_mci <- c(a.plot1[,"mci","soc"], rev(a.plot1[,"mci","int"])) # polygon y-values
yy_mil <- c(a.plot1[,"mil","soc"], rev(a.plot1[,"mil","int"])) # idem
yy_mod <- c(a.plot1[,"mod","soc"], rev(a.plot1[,"mod","int"])) # idem
yy_sev <- c(a.plot1[,"sev","soc"], rev(a.plot1[,"sev","int"])) # idem
yy_dth <- c(a.plot1[,"dth","soc"], rev(a.plot1[,"dth","int"])) # idem
par(mar=c(5, 4, 4, 2)+0.1, xpd=FALSE)
matplot(
  x = v.age_range, 
  y = a.plot1[,,"soc"], 
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
  y = a.plot1[,,"soc"], 
  type = "l",
  lty = 1,
  col = c("green","blue","yellow","orange","black")
)
matlines(
  x = v.age_range, 
  y = a.plot1[,,"int"], 
  type = "l",
  lty = 2,
  col = c("green","blue","yellow","orange","black")
)
legend(x="topright", legend=c("mci","mil","mod","sev","dth"), col=c("green","blue","yellow","orange","black"), lty=1)
legend(x="right", legend=c("soc","int"), col="black", lty=c(1,2))

# plot: mean time in state
par(mar=c(8, 4, 4, 2), xpd=TRUE)
barplot(
  height = m.plot2, 
  horiz = TRUE, 
  xlab = "time (years)", 
  ylab = "strategy", 
  col=c("green","blue","yellow","orange"), 
  space = 0.2, 
  main = "mean time in state"
)
legend(x="bottom", legend=c("mci","mil","mod","sev"), inset=c(0,-0.5), horiz=TRUE, fill=c("green","blue","yellow","orange"))
text(x=c(0,cumsum(m.plot2[1:3,"soc"])), y=1, labels=round(m.plot2[,"soc"],1), pos=4)
text(x=c(0,cumsum(m.plot2[1:3,"int"])), y=2, labels=round(m.plot2[,"int"],1), pos=4)

## table: proportion in state
round(a.plot1[1:10,,"soc"], 3) # state trace standard of care strategy
round(a.plot1[1:10,,"int"], 3) # state trace intervention strategy

## plot: icer
par(mar=c(5, 4, 4, 2)+0.1, xpd=FALSE)
plot(icer, label="all")
## table: icer
print(as.data.frame(t(icer)))

# # optional functions: copy outcomes to clipboard
# 
# ## state trace
# temp.trace_soc <- out_base[["l.out_strategy"]][["soc"]][["m.trace"]]; temp.trace_soc
# write.table(x = temp.trace_soc, file = "clipboard", sep = "\t", row.names = FALSE)
# temp.trace_int <- out_base[["l.out_strategy"]][["int"]][["m.trace"]]; temp.trace_int
# write.table(x = temp.trace_int, file = "clipboard", sep = "\t", row.names = FALSE)
# 
# ## QALYs and costs over time
# temp.out_soc <- out_base[["l.out_strategy"]][["soc"]][["m.out"]]; temp.out_soc
# write.table(x = temp.out_soc, file = "clipboard", sep = "\t", row.names = FALSE)
# temp.out_int <- out_base[["l.out_strategy"]][["int"]][["m.out"]]; temp.out_int
# write.table(x = temp.out_int, file = "clipboard", sep = "\t", row.names = FALSE)

print(m.plot2)


######################################## 5.2.2 DETERMINISTIC SENSITIVITY ANALYSIS ########################################
# to be developed for this version, see release '2.1.0 ISPOR Europe 2023 abstract' for an older version

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

# [to be implemented in updated version]



######################################## 5.2.3 HEADROOM ########################################
# to be developed for this version, see release '2.1.0 ISPOR Europe 2023 abstract' for an older version


