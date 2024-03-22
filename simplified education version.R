
######################################## INFORMATION ########################################

# see readme.md for details



######################################## SOME BASIC R CODING EXPLANATION ########################################

# a vector is an object containing elements
v <- c(11,22,33) # vector containing 3 elements
v # print the vector

# refer to element(s) in a vector using '[]'
v[1] # gives first element of the vector
v[3] # gives third element of the vector
v[c(2,3)] # gives the second and third element of the vector

# replace element(s) in a vector
v[1] <- 9 # replace first element with value 9
v[c(1,2)] <- c(5,6) # replace first and second element with values 5 and 6. Please be aware the length of the elements you replace must be the same (or a multiplicative) of the elements to replace them with
v[] <- 7 # replace all elements with value 7
v[-c(2,3)] <- 5 # replace all elements except the second and third element
v <- 7 # replace the vector with a vector of 1 element; in other words the object 'v' is overwritten
v <- 11:13 # replace the vector with a vector with sequence from 11 to 13

# a matrix is an object with 2 dimensions: rows and columns (opposite to a data frame it can only contain elements of the same type; e.g. numeric or character)
m <- matrix(
  data = c(11,22,33,44,55,66,77,88), # data in the matrix (this in fact is a vector with 8 elements)
  nrow = 4, # number of rows
  ncol = 2, # number of columns
  dimnames = list(
    NULL, # row names
    c("a","b") # column names
  ), 
  byrow = TRUE # fill matrix by row (instead of by column), this determines how the vector (with 8 elements defined under 'data') is put into the matrix
)
m # print the matrix

# refer to element(s) in a matrix
m[1,1] # element in row 1 and column 1
m[3,2] # element in row 3 and column 2
m[3,"b"] # same as previous, although now referring to column name
m[1:4,"b"] # gives elements from row 1 to 4 from column b

# replacing element(s) in a matrix
m[1,1] <- 9 # replace element in row 1 and column 1 with 9
m[,1] <- 9 # replace all elements in column 1 with 9. The single element 9 is recycled to fill all 4 elements in column 1
m[,1] <- c(5,6) # replace all elements in column 1 with 9. The vector with 2 elements (5,6) is recycled to fill all 4 elements in column 1
#m[,1] <- c(5,6,7) # this is not possible as the elements in the matrix to be replaced (4 elements) is not a multiplicative of the elements of the vector to replace them with (3 elements)
m[,] <- 9 # replace all elements with 9
m[,] <- c(1,2) # all elements are replaced with elements 1,2 and are recycled

# array
a <- array(
  data = c(1:24), # data in the array
  dim = c(2,3,4), # length of each dimension, here reflecting 3 elements with length 2, 3 and 4
  dimnames = list( # names of each dimension. a1,a2 are names of the first dimension (can be seen as row names), b1,b2,b3 are names of the second dimension (can be seen as column names), c1,c2,c3,c4 are names of the 3rd dimension
    c("a1","a2"), 
    c("b1","b2","b3"), 
    c("c1","c2","c3","c4")
  )
)
a # print array
a[1,1,1] # element in first location of dimension 1, first location of dimension 2 and first location of dimension 3
a[,,1] # all elements of the third dimension. Because the third dimension only has length of 1, it simplifies the array going from 3 to 2 dimensions
a <- array(data = c(1,2,3,4,11,22,33,44), dim = c(2,2,2))
a[,,1] # all elements of the third dimension
a[1,,] # all elements of the first dimension. Because the third dimension only has length of 1, it simplifies the array going from 3 to 2 dimensions. Notice when printing this it is different from how matrix is printed in its full form

# advanced selection in array
# [to be detailed, perhaps with example of selecting sex- and age-specific death rate]

# list
# [to be detailed]

# loop
# [to be detailed]

# function
# [to be detailed]

# function in function
# [to be detailed]

# debugging
# [to be detailed]

# further read for very basics R: https://rstudio-education.github.io/hopr/basics.html
# further read for advanced R subsetting: https://adv-r.hadley.nz/subsetting.html



######################################## SOME BASIC HEALTH-ECONOMIC MODELING METHODS IN R ########################################

# transition probability table and time dependency operationalized by age-specific mortality
# matrix multiplication
# [to be detailed]



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
m.lifetable_US_2019 <- as.matrix(read.csv(file="life_tables/lifetable_US_2019_ssa.csv", header=TRUE))[,c("male","female")] # import life table and select only men/women and drop age column (make sure age corresponds to row number, i.e., start with age = 1)
m.mortality_rate_US_2019 <- -log(1-(m.lifetable_US_2019)) # convert probability to rate



######################################## 1.2. MODEL INPUTS LIST ########################################

# input parameters (see readme for details)
l.inputs <- list(
  v.names_state = c("mci","mil","mod","sev","dth"), # disease states: mci = mild cognitive impairment; mil = mild dementia; mod = moderate dementia; sev = severe dementia; dth = dead
  v.names_strat = c("soc","int"), # strategies: soc = standard of care; int = intervention
  age_start = 70, 
  sex = "female", 
  p.starting_state_mci = 0.50, 
  n.cycle = 29, 
  p.mci_mil = 0.23, 
  p.mci_mod = 0, 
  p.mci_sev = 0, 
  p.mil_mci = 0.03, 
  p.mil_mod = 0.35, 
  p.mil_sev = 0.04, 
  p.mod_mci = 0, # !!!add this for IPECAD model
  p.mod_mil = 0.03, 
  p.mod_sev = 0.42, 
  p.sev_mci = 0, 
  p.sev_mil = 0, 
  p.sev_mod = 0.02, 
  m.r.mortality = m.mortality_rate_US_2019, 
  hr.mort_mci = 1.82, 
  hr.mort_mil = 2.92, 
  hr.mort_mod = 3.85, 
  hr.mort_sev = 9.52, 
  rr.tx_mci_mil = 0.70, 
  rr.tx_mil_mod = 0.70, 
  u.mci = 0.851 - 0.17, #!!! update these estimates
  u.mil = 0.851 - 0.22, 
  u.mod = 0.851 - 0.36, 
  u.sev = 0.851 - 0.53, 
  c.mci = 6042*1.12 +  460, 
  c.mil = 6042*1.56 +  965 + 0.21*365*0.333, 
  c.mod = 6042*1.93 + 1544 + 0.66*365*0.333, 
  c.sev = 6042*1.93 + 1930, 
  u.Tx_start = -0.14 * (12/52) * 0.035, 
  c.Tx = 26500 + (52/2)*78.35, 
  c.Tx_start = 261.10*4 + 261.10*3*0.215, 
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
    
    # initialize state trace (STEP G3: initialize objects to store strategy outcomes)
    m.trace <- matrix(data = NA, nrow = n.cycle, ncol = n.state, dimnames = list(NULL,v.names_state))
    
    # initialize output table
    m.out <- matrix(data = NA, nrow = n.cycle, ncol = ncol(m.trace)+4, dimnames = list(NULL, c(colnames(m.trace),"ly","qaly","cost","nhb")))
    
    # set starting state distribution (STEP G4: starting state)
    m.trace[1,] <- m.trace1
    
    # markov multiplication (STEP G5: markov multiplication by looping over cycles)
    for(t in 1:(n.cycle-1)) {
      m.trace[t+1,] <- m.trace[t,] %*% a.TP[,,t]
    }
    
    # primary economic outputs (STEP G6: multiply states with utility and cost estimates)
    m.out[,colnames(m.trace)] <- m.trace
    m.out[,"qaly"] <- m.trace %*% c(u.mci, u.mil, u.mod, u.sev, 0) # must match order of states
    m.out[,"cost"] <- m.trace %*% c(c.mci + c.Tx_update, c.mil + c.Tx_update, c.mod, c.sev, 0) # must match order of states
    
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
      m.out[,"qaly"][1] <- m.out[,"qaly"][1] + u.Tx_start
      m.out[,"cost"][1] <- m.out[,"cost"][1] + c.Tx_start
    }
    
    # define vector for discounting QALYs and costs (STEP G8: apply discounting)
    n <- ifelse(test=half_cycle_correction, yes=2, no=1)
    v.discount_EFFECT <- 1 / (( 1 + (discount_EFFECT)) ^ (0 : (n.cycle-n)))
    v.discount_QALY   <- 1 / (( 1 + (discount_QALY))   ^ (0 : (n.cycle-n)))
    v.discount_COST   <- 1 / (( 1 + (discount_COST))   ^ (0 : (n.cycle-n)))
    
    # apply discounting
    m.out[,"ly"] <- m.out[,"ly"] * v.discount_EFFECT
    m.out[,"qaly"] <- m.out[,"qaly"] * v.discount_QALY
    m.out[,"cost"] <- m.out[,"cost"] * v.discount_COST
    
    # calculate net health benefit
    m.out[,"nhb"] <- m.out[,"qaly"] - (m.out[,"cost"] / wtp)
    
    # store strategy-specific output (STEP G9: store outcomes to be wrapped up by the 'run scenario' function)
    return(list(m.out = m.out))
    
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
      
      # death (subset mortality table to obtain age- and sex-specific mortality)
      v.r.dth <- m.r.mortality[age_start:(age_start+n.cycle-1), sex]
      
      # relative risk transitions mci>mil and mil>mod in soc strategy
      rr.mci_mil = 1
      rr.mil_mod = 1
      
      # treatment costs
      c.Tx_update <- 0
      
      # strategy-specific inputs (STEP E: run preparations specific for the intervention strategy)
      if(strat=="int") {
        
        # relative risk transitions mci>mil and mil>mod in int strategy
        rr.mci_mil = rr.tx_mci_mil
        rr.mil_mod = rr.tx_mil_mod
        
        # update transition probabilities treatment effect: during treatment
        p.mci_mil <- 1-exp(-(-log(1-p.mci_mil) * rr.mci_mil)) # convert probability to rate, then multiply with treatment relative risk, then convert to probability
        p.mil_mod <- 1-exp(-(-log(1-p.mil_mod) * rr.mil_mod)) # idem
        
        # treatment costs
        c.Tx_update <- c.Tx
        
      }
      
      # probabilities of remaining in the same state
      p.mci_mci <- 1 - p.mci_mil - p.mci_mod - p.mci_sev
      p.mil_mil <- 1 - p.mil_mci - p.mil_mod - p.mil_sev
      p.mod_mod <- 1 - p.mod_mci - p.mod_mil - p.mod_sev
      p.sev_sev <- 1 - p.sev_mci - p.sev_mil - p.sev_mod
      
      # starting states
      m.trace1 <- matrix(data=0, nrow=1, ncol=n.state, dimnames=list(NULL,v.names_state))
      m.trace1[,"mci"] <- p.starting_state_mci
      m.trace1[,"mil"] <- 1-p.starting_state_mci
      
      # initialize time-dependent TP matrix (STEP G1: prepare transition probability matrix)
      a.TP <- array(data = 0, dim = c(n.state, n.state, n.cycle), dimnames = list(v.names_state,v.names_state,NULL))
      
      # TP matrix state: to death
      a.TP["mci","dth",] <- 1-exp(-(v.r.dth * hr.mort_mci))
      a.TP["mil","dth",] <- 1-exp(-(v.r.dth * hr.mort_mil))
      a.TP["mod","dth",] <- 1-exp(-(v.r.dth * hr.mort_mod))
      a.TP["sev","dth",] <- 1-exp(-(v.r.dth * hr.mort_sev))
      a.TP["dth","dth",] <- 1
      # TP matrix state: from mci
      a.TP["mci","mci",] <- p.mci_mci * (1-a.TP["mci","dth",])
      a.TP["mci","mil",] <- p.mci_mil * (1-a.TP["mci","dth",])
      a.TP["mci","mod",] <- p.mci_mod * (1-a.TP["mci","dth",])
      a.TP["mci","sev",] <- p.mci_sev * (1-a.TP["mci","dth",])
      # TP matrix state: from mild
      a.TP["mil","mci",] <- p.mil_mci * (1-a.TP["mil","dth",])
      a.TP["mil","mil",] <- p.mil_mil * (1-a.TP["mil","dth",])
      a.TP["mil","mod",] <- p.mil_mod * (1-a.TP["mil","dth",])
      a.TP["mil","sev",] <- p.mil_sev * (1-a.TP["mil","dth",])
      # TP matrix state: from moderate
      a.TP["mod","mci",] <- p.mod_mci * (1-a.TP["mod","dth",])
      a.TP["mod","mil",] <- p.mod_mil * (1-a.TP["mod","dth",])
      a.TP["mod","mod",] <- p.mod_mod * (1-a.TP["mod","dth",])
      a.TP["mod","sev",] <- p.mod_sev * (1-a.TP["mod","dth",])
      # TP matrix state: from severe
      a.TP["sev","mci",] <- p.sev_mci * (1-a.TP["sev","dth",])
      a.TP["sev","mil",] <- p.sev_mil * (1-a.TP["sev","dth",])
      a.TP["sev","mod",] <- p.sev_mod * (1-a.TP["sev","dth",])
      a.TP["sev","sev",] <- p.sev_sev * (1-a.TP["sev","dth",])
      
      # check TPs are within 0-1 range
      if(any(a.TP<0 | a.TP>1)) stop("one or more transition probabilities are lower than 0 or higher than 1")
      # check TPs sum to 1 for each cycle (STEP G2: some checks)
      for(i in v.names_state) {
        temp1 <- colSums(a.TP[i,,])
        if(!isTRUE(all.equal(current = temp1, target = rep(1,n.cycle), tolerance = 1e-10))) stop(paste("TPs for",i,"do not add up to 1"))
      }
      
      # list inputs for running each strategy (STEP F: store newly created or updated inputs to be used for the function to run a single strategy)
      l.inputs_strategy <- c(l.inputs, list(
        strat = strat, 
        n.state = n.state, 
        n.strat = n.strat, 
        n.cycle = n.cycle, 
        c.Tx_update = c.Tx_update, 
        a.TP = a.TP, 
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
    
    # return result
    if(detailed) {
      return(list(
        df.out = df.out,
        l.out = l.out
      ))
    }
    if(!detailed) return(df.out)
  }
  )
}



######################################## 3. MODEL CALIBRATION ########################################
# n/a



######################################## 4. VALIDATION ########################################
# n/a



######################################## 5. ANALYSIS ########################################

######################################## 5.1. DETERMINISTIC ########################################

# run scenario and results
l.out <- f.run_scenario(l.inputs = l.inputs, detailed = TRUE)
print(l.out$df.out)

# standard tables/plots
if(F) {
  
  # all outcomes
  # str(l.out) # show structure of output
  # l.out[["df.out"]]
  # l.out # all detailed output
  # l.out[["l.out"]] # strategy details
  # l.out[["l.out"]][["soc"]] # strategy 'soc' details
  # l.out[["l.out"]][["int"]] # strategy 'int' details
  # l.out[["l.out"]][["soc"]][["a.TP"]] # strategy 'soc' transition probability matrix
  # l.out[["l.out"]][["soc"]][["m.trace"]] # strategy 'soc' state trace
  # l.out[["l.out"]][["soc"]][["m.out"]] # strategy 'soc' outcomes per cycle
  # l.out[["l.out"]][["int"]][["a.TP"]] # idem for 'int'
  # l.out[["l.out"]][["int"]][["m.trace"]]
  # l.out[["l.out"]][["int"]][["m.out"]]
  
  # table: summary
  m.result <- matrix(
    data = NA, 
    nrow = ncol(l.out$l.out$soc$m.out), 
    ncol = 3, 
    dimnames = list(
      colnames(l.out$l.out$soc$m.out),
      c("soc","int","dif")
    )
  )
  
  # fill matrix
  for(i in rownames(m.result)) {
    m.result[i,"soc"] <- sum(l.out$l.out$soc$m.out[,i])
    m.result[i,"int"] <- sum(l.out$l.out$int$m.out[,i])
  }
  m.result[,"dif"] <- m.result[,"int"] - m.result[,"soc"]
  
  # plot: time in state
  plot.timestate_data <- m.result[c("mci","mil","mod","sev"),c("int","soc")]
  
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
  l.out$l.out$soc$m.out[,c("mci","mil","mod","sev","dth")]
  l.out$l.out$int$m.out[,c("mci","mil","mod","sev","dth")]
  
  # plot: state trace
  v.age_range <- c(l.inputs[["age_start"]]:(l.inputs[["age_start"]]+l.inputs[["n.cycle"]]-2)) # store age range
  xx <- c(v.age_range, rev(v.age_range)) # prepare polygon x-values
  yy_mci <- c(l.out$l.out$soc$m.out[,"mci"], rev(l.out$l.out$int$m.out[,"mci"])) # polygon y-values
  yy_mil <- c(l.out$l.out$soc$m.out[,"mil"], rev(l.out$l.out$int$m.out[,"mil"])) # idem
  yy_mod <- c(l.out$l.out$soc$m.out[,"mod"], rev(l.out$l.out$int$m.out[,"mod"])) # idem
  yy_sev <- c(l.out$l.out$soc$m.out[,"sev"], rev(l.out$l.out$int$m.out[,"sev"])) # idem
  yy_dth <- c(l.out$l.out$soc$m.out[,"dth"], rev(l.out$l.out$int$m.out[,"dth"])) # idem
  
  par(mar=c(5, 4, 4, 1)+0.1, xpd=FALSE)
  matplot(
    x = v.age_range, 
    y = cbind(l.out$l.out$soc$m.out[,c("mci","mil","mod","sev","dth")],l.out$l.out$int$m.out[,c("mci","mil","mod","sev","dth")]), 
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
    y = l.out$l.out$soc$m.out[,c("mci","mil","mod","sev","dth")], 
    type = "l",
    lty = 1,
    col = c("green","yellow","orange","red","black")
  )
  matlines(
    x = v.age_range, 
    y = l.out$l.out$int$m.out[,c("mci","mil","mod","sev","dth")], 
    type = "l",
    lty = 2,
    col = c("green","yellow","orange","red","black")
  )
  legend(x="topright", legend=c("mci","mil","mod","sev","dth"), col=c("green","yellow","orange","red","black"), lty=1, bg="white")
  legend(x="right", legend=c("soc","int"), col="black", lty=c(1,2), bg="white")
  
  # table: incremental cost-effectiveness ratio
  icer <- calculate_icers(cost = l.out[["df.out"]][,"COST"], effect = l.out[["df.out"]][,"QALY"], strategies = l.out[["df.out"]][,"strategy"])
  as.data.frame(icer)
  
  # plot: incremental cost-effectiveness plane
  par(mar=c(5, 4, 4, 1)+0.1, xpd=FALSE)
  print(plot(icer, label="all"))
  
}



######################################## 5.1. DETERMINISTIC SENSITIVITY ANALYSIS ########################################

if(T) {
  
  # list parameters for DSA
  dsa_pars <- c(
    "age_start",
    "p.mci_mil",
    "p.mil_mod",
    "hr.mort_mil",
    "rr.tx_mci_mil",
    "rr.tx_mil_mod",
    "u.mci",
    "c.mci",
    "c.Tx"
  )
  
  # list minimum values
  dsa_min <- c(
    l.inputs[["age_start"]] - 10,
    l.inputs[["p.mci_mil"]] - 0.05,
    l.inputs[["p.mil_mod"]] - 0.05,
    l.inputs[["hr.mort_mil"]] - 0.3,
    l.inputs[["rr.tx_mci_mil"]] - 0.10,
    l.inputs[["rr.tx_mil_mod"]] - 0.10,
    l.inputs[["u.mci"]] - 0.1,
    l.inputs[["c.mci"]] - 0.1,
    l.inputs[["c.Tx"]] - 1000
  )
  
  # list maximum values
  dsa_max <- c(
    l.inputs[["age_start"]] + 10,
    l.inputs[["p.mci_mil"]] + 0.05,
    l.inputs[["p.mil_mod"]] + 0.05,
    l.inputs[["hr.mort_mil"]] + 0.3,
    l.inputs[["rr.tx_mci_mil"]] + 0.10,
    l.inputs[["rr.tx_mil_mod"]] + 0.10,
    l.inputs[["u.mci"]] + 0.1,
    l.inputs[["c.mci"]] + 0.1,
    l.inputs[["c.Tx"]] + 1000
  )
  
  # define sensitivity analysis range
  df.owsa_params_range <- data.frame(
    pars = dsa_pars,
    min = dsa_min,
    max = dsa_max
  )
  df.owsa_params_range
  
  # run DSA
  owsa_det <- run_owsa_det(
    params_range = df.owsa_params_range,
    params_basecase = l.inputs,
    nsamp = 10, # n equally-spaced samples
    FUN = f.run_scenario, # make sure the function only returns the aggregated outcomes (not a list of outcomes)
    outcomes = NULL, # NULL for all outcomes (which should then be selected when plotting a tornado graph)
    strategies = NULL, # all strategies
    progress = TRUE
  )
  #str(owsa_det)
  #owsa_det
  
  # Plot outcome of each strategy over each parameter range
  plot(owsa_det$owsa_NHB, n_x_ticks = 3)
  
}