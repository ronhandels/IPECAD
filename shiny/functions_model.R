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
