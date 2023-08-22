


######################################## INFORMATION ########################################

# see readme for details



######################################## TECHNICAL PREPARATION ########################################

cat("\014") # clear console
rm(list = ls()) # clear environment
# install.packages("dampack") # remove # and run once to install package
library(dampack) # load package
setwd("~/GitHub/IPECAD")
#setwd("D:/surfdrive/PhD/Projects/IPECAD/open source model/2.0/_github/IPECAD/") # set working directory; change to your own directory
#setwd("C:/users/Ron/surfdrive/PhD/Projects/IPECAD/open source model/2.0/_github/IPECAD/") # set working directory; change to your own directory


######################################## 1. INPUTS ########################################

######################################## 1.1. ESTIMATED ########################################

# U.S. general population life table
m.lifetable_US <- as.matrix(read.csv(file="life_tables/lifetable_US.csv", header=TRUE))[,c("male","female")] # import life table and select only men/women and drop age column (make sure age corresponds to row number, i.e., start with age = 1)
m.mortality_rate_US <- -log(1-(m.lifetable_US)) # convert probability to rate



######################################## 1.2. MODEL INPUTS LIST ########################################

l.inputs <- list(
  v.names_state = c("mcion","mciof","milon","milof","mod","sev","dth"), # disease states: mci = mild cognitive impairment; mil = mild dementia; mod = moderate dementia; sev = severe dementia; dth = dead
  v.names_strat = c("soc","int"), # strategies: soc = standard of care strategy; int = intervention strategy
  age_start = 70, # age of starting population
  age_end = 100, # age up to which model is run (reflecting time horizon)
  sex = "female", # sex of starting population
  p.mci_mil = 0.206, # p.x_x: transition probability between states [Wimo, 2020: https://doi.org/10.3233/jad-191055]
  p.mci_mod = 0, # idem
  p.mci_sev = 0, # idem
  p.mil_mci = 0, # idem
  p.mil_mod = 0.293, # idem
  p.mil_sev = 0.001, # idem
  p.mod_mil = 0.087, # idem
  p.mod_sev = 0.109, # idem
  p.sev_mil = 0.000, # idem
  p.sev_mod = 0.196, # idem
  m.r.mortality = m.mortality_rate_US, # general population mortality rate
  hr.mort_mci = 1, # hr.mort_x: hazard ratio mortality by disease state (hazard ratio by dementia state compared to very mild dementia [Wimo, 2020: https://doi.org/10.3233/jad-191055] multiplied with HR or very mild compared to no dementia [Andersen, 2010: https://doi.org/10.1159/000265553])
  hr.mort_verymilddem = 1.82, # [Andersen, 2010: https://doi.org/10.1159/000265553]
  hr.mort_mil = 1.318, # [Wimo, 2020: https://doi.org/10.3233/jad-191055]
  hr.mort_mod = 2.419, # idem
  hr.mort_sev = 4.267, # idem
  rr.tx_mci_mil = 0.75, # treatment effect expressed as hazard ratio on transition rate
  rr.tx_mci_mod = 1, # idem
  rr.tx_mci_sev = 1, # idem
  rr.tx_mil_mod = 0.75, # assumed same effect as for MCI to dementia
  tx_waning = 0.05, # assumed annual waning of treatment
  p.discontinuation1 = 0.1, # discontinuation at year 1
  p.discontinuation_x = 0.1, # annual proportion discontinuation after year 1
  tx_duration = 7, # maximum treatment duration
  p.starting_state_mci = 1, # proportion starting in disease state MCI, remaining from 1 will start in 'mil' (all will start as 'of' in 'soc' and 'on' in 'int')
  u.mci = 0.73, # u.x: utility in state [https://doi.org/10.1016/j.jalz.2019.05.004]
  u.mil = 0.69, # idem
  u.mod = 0.53, # idem
  u.sev = 0.38, # idem
  c.mci = (1254 +  222) * 12 * (1-0    ) + (1254 + 8762) * 12 * 0, # c.x: costs in state, build up as montly costs in patient health and social care by care setting (community/residential) multiplied by 12 (annual costs) [Tahami, 2023: https://doi.org/10.1007/s40120-023-00460-1 table 2] and multiplied by proportion in setting [Tahami, 2023: https://doi.org/10.1007/s40120-023-00460-1 table 1]
  c.mil = (1471 +  410) * 12 * (1-0.038) + (1471 + 8762) * 12 * 0.038, # idem
  c.mod = (1958 +  653) * 12 * (1-0.110) + (1958 + 8762) * 12 * 0.110, # idem
  c.sev = (2250 + 1095) * 12 * (1-0.259) + (2250 + 8762) * 12 * 0.259, # idem
  c.Tx = 10000, # treatment costs
  c.Tx_diagnostics1 = 2000, # costs diagnostics cycle 1 (not half-cycle corrected)
  discount_QALY = 0.035, # discount rate
  discount_COST = 0.035, # # discount rate
  wtp = 40000 # willingness to pay
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
  # G7: apply half-cycle correction
  # G8: apply discounting
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
    a.TP["mod","dth",] <- 1-exp(-(v.r.dth * hr.mort_mod * hr.mort_verymilddem))
    a.TP["sev","dth",] <- 1-exp(-(v.r.dth * hr.mort_sev * hr.mort_verymilddem))
    a.TP["dth","dth",] <- 1
    
    # TP matrix state: from mci-on
    a.TP["mcion","mcion",] <- v.p.mcion_mci * (1-v.p.discontinuation) * (1-a.TP["mcion","dth",])
    a.TP["mcion","mciof",] <- v.p.mcion_mci *    v.p.discontinuation  * (1-a.TP["mcion","dth",])
    a.TP["mcion","milon",] <- v.p.mcion_mil * (1-v.p.discontinuation) * (1-a.TP["mcion","dth",])
    a.TP["mcion","milof",] <- v.p.mcion_mil *    v.p.discontinuation  * (1-a.TP["mcion","dth",])
    a.TP["mcion","mod",]   <- v.p.mcion_mod                           * (1-a.TP["mcion","dth",])
    a.TP["mcion","sev",]   <- v.p.mcion_sev                           * (1-a.TP["mcion","dth",])

    # TP matrix state: from mci-off
    a.TP["mciof","mciof",] <- v.p.mci_mci                           * (1-a.TP["mciof","dth",])
    a.TP["mciof","milof",] <- v.p.mci_mil                           * (1-a.TP["mciof","dth",])
    a.TP["mciof","mod",]   <- v.p.mci_mod                           * (1-a.TP["mcion","dth",])
    a.TP["mciof","sev",]   <- v.p.mci_sev                           * (1-a.TP["mcion","dth",])
    
    # TP matrix state: from mild-on
    a.TP["milon","mcion",] <- v.p.milon_mci * (1-v.p.discontinuation) * (1-a.TP["milon","dth",])
    a.TP["milon","mciof",] <- v.p.milon_mci *    v.p.discontinuation  * (1-a.TP["milon","dth",])
    a.TP["milon","milon",] <- v.p.milon_mil * (1-v.p.discontinuation) * (1-a.TP["milon","dth",])
    a.TP["milon","milof",] <- v.p.milon_mil *    v.p.discontinuation  * (1-a.TP["milon","dth",])
    a.TP["milon","mod",]   <- v.p.milon_mod                           * (1-a.TP["milon","dth",])
    a.TP["milon","sev",]   <- v.p.mil_sev                             * (1-a.TP["milon","dth",])
    
    # TP matrix state: from mild-off
    a.TP["milof","mciof",] <- v.p.mil_mci                           * (1-a.TP["milof","dth",])
    a.TP["milof","milof",] <- v.p.mil_mil                           * (1-a.TP["milof","dth",])
    a.TP["milof","mod",]   <- v.p.mil_mod                           * (1-a.TP["milof","dth",])
    a.TP["milof","sev",]   <- v.p.mil_sev                           * (1-a.TP["milof","dth",])

    # TP matrix state: from moderate
    a.TP["mod","milof",] <- v.p.mod_mil * (1-a.TP["mod","dth",])
    a.TP["mod","mod",]   <- v.p.mod_mod * (1-a.TP["mod","dth",])
    a.TP["mod","sev",]   <- v.p.mod_sev * (1-a.TP["mod","dth",])
    
    # TP matrix state: from severe
    a.TP["sev","milof",] <- v.p.sev_mil * (1-a.TP["sev","dth",])
    a.TP["sev","mod",]   <- v.p.sev_mod * (1-a.TP["sev","dth",])
    a.TP["sev","sev",]   <- v.p.sev_sev * (1-a.TP["sev","dth",])
    
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
    m.out[,"ly"]   <- m.trace %*% c(1           , 1    , 1           , 1    , 1    , 1    , 0) # must match order of states
    m.out[,"qaly"] <- m.trace %*% c(u.mci       , u.mci, u.mil       , u.mil, u.mod, u.sev, 0) # must match order of states
    m.out[,"cost"] <- m.trace %*% c(c.mci + c.Tx, c.mci, c.mil + c.Tx, c.mil, c.mod, c.sev, 0) # must match order of states
    
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
      
      # treatment costs
      c.Tx <- 0
      
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

if(F) {
  
  # set number of PSA replications
  n.psa <- 1000
  
  # list parameters for PSA
  psa_params <- c(
    "age_start",
    "hr.mort_verymilddem",
    "hr.mort_mil",
    "hr.mort_mod",
    "hr.mort_sev",
    "p.mci_mil",
    "p.mil_mod",
    "p.mil_sev",
    "p.mod_mil",
    "p.mod_sev",
    "p.sev_mil",
    "p.sev_mod",
    "rr.tx_mci_mil",
    "rr.tx_mil_mod",
    "tx_duration",
    "c.Tx"
  )
  
  # list distributions
  psa_dists <- c(
    "bootstrap",
    "log-normal",
    "log-normal",
    "log-normal",
    "log-normal",
    "truncated-normal",
    "bootstrap",
    "bootstrap",
    "bootstrap",
    "bootstrap",
    "bootstrap",
    "bootstrap",
    "truncated-normal",
    "truncated-normal",
    "bootstrap",
    "bootstrap"
  )
  
  # parameterization types
  psa_parameterization_types <- c(
    "value, weight",
    "mean, sd",
    "mean, sd",
    "mean, sd",
    "mean, sd",
    "mean, sd, ll, ul",
    "value, weight",
    "value, weight",
    "value, weight",
    "value, weight",
    "value, weight",
    "value, weight",
    "mean, sd, ll, ul",
    "mean, sd, ll, ul",
    "value, weight",
    "value, weight"
  )
  
  # sample transition probabilities based on oprobit
  ## function to produce transition probability matrix using ordered probit regression coefficients
  f.oprobit_TP <- function(b_mod, b_sev, cut1, cut2) {
    p.mil_mil <- pnorm(cut1 - 0) - pnorm(-999 - 0) # mil>mil
    p.mil_mod <- pnorm(cut2 - 0) - pnorm(cut1 - 0) # mil>mod
    p.mil_sev <- pnorm(999 - 0) - pnorm(cut2 - 0) # mil>mod
    p.mod_mil <- pnorm(cut1 - b_mod) - pnorm(-999 - b_mod) # mod>mil
    p.mod_mod <- pnorm(cut2 - b_mod) - pnorm(cut1 - b_mod) # mod>mod
    p.mod_sev <- pnorm(999 - b_mod) - pnorm(cut2 - b_mod) # mod>mod
    p.sev_mil <- pnorm(cut1 - b_sev) - pnorm(-999 - b_sev) # sev>mil
    p.sev_mod <- pnorm(cut2 - b_sev) - pnorm(cut1 - b_sev) # sev>mod
    p.sev_sev <- pnorm(9999 - b_sev) - pnorm(cut2 - b_sev) # sev>mod
    return(cbind(p.mil_mil, p.mil_mod, p.mil_sev, p.mod_mil, p.mod_mod, p.mod_sev, p.sev_mil, p.sev_mod, p.sev_sev))
  }
  
  # ## general base case transition probability matrix
  # round(f.oprobit_TP(b_mod=1.8984, b_sev=3.9837, cut1=0.542, cut2=3.129), 3)
  
  ## sample ordered probit parameter values (using 'rnorm') and convert ordered probit parameters to transition probabilities (using 'f.oprobit_TP')
  m.TP_psa <- f.oprobit_TP(
    b_mod = rnorm(n=n.psa, mean=1.8984, sd=(1.949-1.8479)/(1.96*4)), # 95% confidence interval upper bound minus lower bound, then divide by 4 times standard deviation
    b_sev = rnorm(n=n.psa, mean=3.9837, sd=(4.233-3.7343)/(1.96*4)),
    cut1 = rnorm(n=n.psa, mean=0.542, sd=(0.574-0.510)/(1.96*4)),
    cut2 = rnorm(n=n.psa, mean=3.129, sd=(3.206-3.052)/(1.96*4))
  ) # m.TP_psa # hist(m.TP_psa[,2])
  
  # distribution parameters (seems that weights must add up to 1 to be sampled from their explicit values rather than be seen as a distribution to be sampled from, see help)
  psa_dists_params <- list(
    data.frame(value=c(60,65,70,75,80), weight=c(0.1,0.2,0.4,0.2,0.1)), # data.frame(value=c(60,65,70,75,80), weight=c(0.5,0.75,1,0.75,0.5))
    c(l.inputs[["hr.mort_verymilddem"]], (2.14  - 1.55 )/(qnorm(0.975)*2)), # 'qnorm(0.975)*2' represents 95%CI interval (i.e., 4 times standard deviation)
    c(l.inputs[["hr.mort_mil"]]        , (1.507 - 1.153)/(qnorm(0.975)*2)),
    c(l.inputs[["hr.mort_mod"]]        , (2.757 - 2.122)/(qnorm(0.975)*2)),
    c(l.inputs[["hr.mort_sev"]]        , (5.043 - 3.610)/(qnorm(0.975)*2)),
    c(l.inputs[["p.mci_mil"]], 0.05, 0, 1),
    data.frame(value=m.TP_psa[,"p.mil_mod"], weight=rep(1,n.psa)),
    data.frame(value=m.TP_psa[,"p.mil_sev"], weight=rep(1,n.psa)),
    data.frame(value=m.TP_psa[,"p.mod_mil"], weight=rep(1,n.psa)),
    data.frame(value=m.TP_psa[,"p.mod_sev"], weight=rep(1,n.psa)),
    data.frame(value=m.TP_psa[,"p.sev_mil"], weight=rep(1,n.psa)),
    data.frame(value=m.TP_psa[,"p.sev_mod"], weight=rep(1,n.psa)),
    c(l.inputs[["rr.tx_mci_mil"]], 0.10, 0, 1),
    c(l.inputs[["rr.tx_mci_mil"]], 0.10, 0, 1),
    data.frame(value=c(7,10,15), weight=c(0.2,0.6,0.2)),
    data.frame(value=c(1000,5000,10000,20000,25000), weight=c(0.1,0.2,0.4,0.2,0.1))
  )
  
  # sample parameter values based on underlying distributions independently
  out_psa_samp <- gen_psa_samp(
    params = psa_params,
    dists = psa_dists,
    parameterization_types = psa_parameterization_types,
    dists_params = psa_dists_params,
    n = n.psa
  ) # ; out_psa_samp; hist(out_psa_samp$c.Tx)
  
  # run model using sampled input parameters
  out_psa <- run_psa(
    psa_samp = out_psa_samp,
    params_basecase = l.inputs,
    FUN = f.run_scenario,
    outcomes = c("QALY", "COST"),
    strategies = c("soc", "int"),
    progress = TRUE
  )
  
  # make PSA object
  obj_psa <- make_psa_obj(
    cost = out_psa$COST$other_outcome,
    effect = out_psa$QALY$other_outcome,
    parameters = out_psa$COST$parameters,
    strategies = out_psa$COST$strategies,
    currency = "$"
  )
  
  # ICE-plane
  plot(obj_psa)
  
  # CEAC
  obj_ceac <- ceac(wtp = seq(from=0, to=150000, by=10000), psa = obj_psa)
  plot(obj_ceac, frontier = TRUE, points = TRUE)
  
  # OWSA
  owsa_psa <- owsa(sa_obj=obj_psa, outcome="eff") # create owsa on all parameters (see dampack vignette for details)
  owsa_tornado(owsa_psa, n_y_ticks = 2) # plot tornado graph
  owsa_psa <- owsa(sa_obj=obj_psa, outcome="cost") # create owsa on all parameters (see dampack vignette for details)
  owsa_tornado(owsa_psa, n_y_ticks = 2) # plot tornado graph
  
}



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
  mci=out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"mcion"] + out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"mciof"], 
  mil=out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"milon"] + out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"milof"], 
  mod=out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"mod"], 
  sev=out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"sev"], 
  dth=out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"dth"]
)
m.plot1_int <- cbind(
  mci=out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"mcion"] + out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"mciof"], 
  mil=out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"milon"] + out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"milof"], 
  mod=out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"mod"], 
  sev=out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"sev"], 
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
  col = c("green","yellow","orange","red","black")
)
matlines(
  x = v.age_range, 
  y = a.plot1[,,"int"], 
  type = "l",
  lty = 2,
  col = c("green","yellow","orange","red","black")
)
legend(x="topright", legend=c("mci","mil","mod","sev","dth"), col=c("green","yellow","orange","red","black"), lty=1)
legend(x="right", legend=c("soc","int"), col="black", lty=c(1,2))

# plot: mean time in state
par(mar=c(8, 4, 4, 2), xpd=TRUE)
barplot(
  height = m.plot2, 
  horiz = TRUE, 
  xlab = "time (years)", 
  ylab = "strategy", 
  col=c("green","yellow","orange","red"), 
  space = 0.2, 
  main = "mean time in state"
)
legend(x="bottom", legend=c("mci","mil","mod","sev"), inset=c(0,-0.5), horiz=TRUE, fill=c("green","yellow","orange","red"))
text(x=c(0,cumsum(m.plot2[1:3,"soc"])), y=1, labels=round(m.plot2[,"soc"],1), pos=4)
text(x=c(0,cumsum(m.plot2[1:3,"int"])), y=2, labels=round(m.plot2[,"int"],1), pos=4)

## table: proportion in state
round(a.plot1[1:10,,"soc"], 2) # state trace standard of care strategy
round(a.plot1[1:10,,"int"], 2) # state trace intervention strategy

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



######################################## 5.2.2 DETERMINISTIC SENSITIVITY ANALYSIS ########################################

if(F) {
  
  # list parameters for DSA
  dsa_pars <- c(
    "age_start",
    "hr.mort_verymilddem",
    "hr.mort_mil",
    "hr.mort_mod",
    "hr.mort_sev",
    "p.mci_mil",
    "p.mil_mod",
    "p.mod_mil",
    "p.mod_sev",
    "p.sev_mod",
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
    quantile(x=m.TP_psa[,"p.mil_mod"], probs=0.025),
    quantile(x=m.TP_psa[,"p.mod_mil"], probs=0.025),
    quantile(x=m.TP_psa[,"p.mod_sev"], probs=0.025),
    quantile(x=m.TP_psa[,"p.sev_mod"], probs=0.025),
    l.inputs[["rr.tx_mci_mil"]]^2,
    l.inputs[["rr.tx_mil_mod"]]^2,
    l.inputs[["tx_duration"]]/2,
    l.inputs[["c.Tx"]]/2
  )
  
  # list maximum values
  dsa_max <- c(
    80,
    2.14,
    1.507,
    2.757,
    5.043,
    l.inputs[["p.mci_mil"]]*2,
    quantile(x=m.TP_psa[,"p.mil_mod"], probs=0.975),
    quantile(x=m.TP_psa[,"p.mod_mil"], probs=0.975),
    quantile(x=m.TP_psa[,"p.mod_sev"], probs=0.975),
    quantile(x=m.TP_psa[,"p.sev_mod"], probs=0.975),
    l.inputs[["rr.tx_mci_mil"]]^0.5,
    l.inputs[["rr.tx_mil_mod"]]^0.5,
    l.inputs[["tx_duration"]]*2,
    l.inputs[["c.Tx"]]*2
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
    nsamp = 20, # n equally-spaced samples
    FUN = f.run_scenario, # make sure the function only returns the aggregated outcomes (not a list of outcomes)
    outcomes = NULL, # NULL for all outcomes (which should then be selected when plotting a tornado graph)
    strategies = NULL, # all strategies
    progress = TRUE
  )

  # OWSA: plot tornado graph of DSA outcome 'NHB' (only 1 outcome is possible; possibly can't work with negative (NHB) data?)
  owsa_tornado(out_owsa_det$owsa_QALY)
  
  # OWSA: plot tornado graph of DSA outcome 'iNHB'
  out_owsa_det$owsa_NHB
  out_owsa_det_NHB_soc <- subset(x=out_owsa_det$owsa_NHB, subset=strategy=="st_soc")
  out_owsa_det_NHB_int <- subset(x=out_owsa_det$owsa_NHB, subset=strategy=="st_int")
  if(!all(out_owsa_det_NHB_soc[,c("parameter","param_val")]==out_owsa_det_NHB_int[,c("parameter","param_val")])) stop("error")
  out_owsa_det_iNHB <- out_owsa_det_NHB_int
  out_owsa_det_iNHB[,"outcome_val"] <- NA
  out_owsa_det_iNHB[,"outcome_val"] <- out_owsa_det_NHB_int[,"outcome_val"] - out_owsa_det_NHB_soc[,"outcome_val"] + 100 # added 100 because negative values produces an error
  owsa_tornado(out_owsa_det_iNHB, n_y_ticks=4)
  
}



######################################## 5.2.3 HEADROOM ########################################

if(F) {
  
  # import all life tables
  a.lifetable <- array(data=NA, dim=c(100,2,6), dimnames=list(NULL,c("male","female"),c("ES","NL","PL","SE","UK","US")))
  a.lifetable[,,"ES"] <- as.matrix(read.csv(file="life_tables/lifetable_ES.csv", header=TRUE))[,c("male","female")]
  a.lifetable[,,"NL"] <- as.matrix(read.csv(file="life_tables/lifetable_NL.csv", header=TRUE))[,c("male","female")]
  a.lifetable[,,"PL"] <- as.matrix(read.csv(file="life_tables/lifetable_PL.csv", header=TRUE))[,c("male","female")]
  a.lifetable[,,"SE"] <- as.matrix(read.csv(file="life_tables/lifetable_SE.csv", header=TRUE))[,c("male","female")]
  a.lifetable[,,"UK"] <- as.matrix(read.csv(file="life_tables/lifetable_UK.csv", header=TRUE))[,c("male","female")]
  a.lifetable[,,"US"] <- as.matrix(read.csv(file="life_tables/lifetable_US.csv", header=TRUE))[,c("male","female")]
  matplot(x=a.lifetable[1:99,"male",], type="l", col=rainbow(6))
  legend(x="topleft", legend=c("ES","NL","PL","SE","UK","US"), col=rainbow(6), lty=c(1:6))
  
  # costs (Euro) in EU regions (source: Jonsson, 2023 https://doi.org/10.1007/s40273-022-01212-z table 4 consumer price index is 2021)
  m.c_region <- matrix(data=c(
    NA, NA, NA, NA, NA, 
    19909, 7616, 20876, 20420, 31984, 
    34223, 9670, 37540, 40953, 47934, 
    61958, 11236, 58198, 61906, 56104
  ), 
  nrow=5, 
  ncol=4, 
  dimnames = list(c("britishisles","eastbaltics","north","south","west"),c("mci","mil","mod","sev"))
  ); m.c_region
  m.c_region[,"mci"] <- m.c_region[,"mil"] * (17712/26380.51) # add costs for MCI (assumed ratio between m.mil and m.mci from base case inputs)
  
  # quality of life
  l.inputs_headroom <- l.inputs
  l.inputs_headroom[["u.mci"]] <- mean(c(0.7 ,0.75)) # Landeiro review 2020 studies Jonsson 2006 and Hessman 2016 (as these are EU mixed setting including MCI and dementia stages)
  l.inputs_headroom[["u.mil"]] <- mean(c(0.65,0.61))
  l.inputs_headroom[["u.mod"]] <- mean(c(0.51,0.41))
  l.inputs_headroom[["u.sev"]] <- mean(c(0.40,0.21))
  
  # costs
  l.inputs_UK <- l.inputs_headroom
  l.inputs_UK[["c.mci"]] <- m.c_region["britishisles","mci"]
  l.inputs_UK[["c.mil"]] <- m.c_region["britishisles","mil"]
  l.inputs_UK[["c.mod"]] <- m.c_region["britishisles","mod"]
  l.inputs_UK[["c.sev"]] <- m.c_region["britishisles","sev"]
  l.inputs_UK[["m.r.mortality"]] <- -log(1-(a.lifetable[,,"UK"]))
  l.inputs_E <- l.inputs_headroom
  l.inputs_E[["c.mci"]] <- m.c_region["eastbaltics","mci"]
  l.inputs_E[["c.mil"]] <- m.c_region["eastbaltics","mil"]
  l.inputs_E[["c.mod"]] <- m.c_region["eastbaltics","mod"]
  l.inputs_E[["c.sev"]] <- m.c_region["eastbaltics","sev"]
  l.inputs_E[["m.r.mortality"]]  <- -log(1-(a.lifetable[,,"PL"]))
  l.inputs_N <- l.inputs_headroom
  l.inputs_N[["c.mci"]] <- m.c_region["north","mci"]
  l.inputs_N[["c.mil"]] <- m.c_region["north","mil"]
  l.inputs_N[["c.mod"]] <- m.c_region["north","mod"]
  l.inputs_N[["c.sev"]] <- m.c_region["north","sev"]
  l.inputs_N[["m.r.mortality"]]  <- -log(1-(a.lifetable[,,"SE"]))
  l.inputs_S <- l.inputs_headroom
  l.inputs_S[["c.mci"]] <- m.c_region["south","mci"]
  l.inputs_S[["c.mil"]] <- m.c_region["south","mil"]
  l.inputs_S[["c.mod"]] <- m.c_region["south","mod"]
  l.inputs_S[["c.sev"]] <- m.c_region["south","sev"]
  l.inputs_S[["m.r.mortality"]]  <- -log(1-(a.lifetable[,,"ES"]))
  l.inputs_W <- l.inputs_headroom
  l.inputs_W[["c.mci"]] <- m.c_region["west","mci"]
  l.inputs_W[["c.mil"]] <- m.c_region["west","mil"]
  l.inputs_W[["c.mod"]] <- m.c_region["west","mod"]
  l.inputs_W[["c.sev"]] <- m.c_region["west","sev"]
  l.inputs_W[["m.r.mortality"]]  <- -log(1-(a.lifetable[,,"NL"]))
  l.inputs_US <- l.inputs

  # headroom function
  f.headroom <- function(x, l.inputs, parameter) {
    l.inputs[[parameter]] <- x
    out <- f.run_scenario(l.inputs=l.inputs, detailed=FALSE)
    iNHB <- out[2,"NHB"] - out[1,"NHB"]
    return(iNHB)
  }
  
  # run headroom
  headroom_UK <- optimize(f=function(x) abs(f.headroom(x, l.inputs=l.inputs_UK, parameter="c.Tx")), interval=c(1,100000))[["minimum"]]; headroom_UK
  headroom_E <- optimize(f=function(x) abs(f.headroom(x, l.inputs=l.inputs_E, parameter="c.Tx")), interval=c(1,100000))[["minimum"]]; headroom_E
  headroom_N <- optimize(f=function(x) abs(f.headroom(x, l.inputs=l.inputs_N, parameter="c.Tx")), interval=c(1,100000))[["minimum"]]; headroom_N
  headroom_S <- optimize(f=function(x) abs(f.headroom(x, l.inputs=l.inputs_S, parameter="c.Tx")), interval=c(1,100000))[["minimum"]]; headroom_S
  headroom_W <- optimize(f=function(x) abs(f.headroom(x, l.inputs=l.inputs_W, parameter="c.Tx")), interval=c(1,100000))[["minimum"]]; headroom_W
  headroom_US <- optimize(f=function(x) abs(f.headroom(x, l.inputs=l.inputs_US, parameter="c.Tx")), interval=c(1,100000))[["minimum"]]; headroom_US
  
  # result
  round(c(headroom_UK, headroom_E, headroom_N, headroom_S, headroom_W, headroom_US)) # difference with ISPOR 2023 abstract due to patch that changed the order of discounting and half-cycle correction
  
}



