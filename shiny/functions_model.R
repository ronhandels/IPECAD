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
    
    # check TPs are within 0-1 range
    if(any(a.TP<0 | a.TP>1)) stop("one or more transition probabilities are lower than 0 or higher than 1")
    # check TPs sum to 1 for each cycle (STEP G2: some checks)
    for(i in v.names_state) {
      temp1 <- colSums(a.TP[i,,])
      if(!isTRUE(all.equal(current = temp1, target = rep(1,n.cycle), tolerance = 1e-10))) stop(paste("TPs for",i,"do not add up to 1"))
    }
    
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
      m.out <- m.out[-n.cycle,] # remove the last cycle
    }
    
    # add additional inputs to cycle 1
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
    
    # some validity checks on input estimates
    # [to be developed]
    
    # store counters (STEP B: prepare and initialize objects to store scenario and strategy outcomes)
    n.state <- length(v.names_state) # number of states
    n.strat <- length(v.names_strat) # number of strategies
    
    # initialize output data frame (create an empty data frame to store outcomes of a scenario)
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




