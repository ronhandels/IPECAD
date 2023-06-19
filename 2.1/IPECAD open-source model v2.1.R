

# to change: 
# . additional costs after half-cycle correction (or double the costs)
# . update utilities and cost estimates
# . produce IPECAD cross-comparison outcomes (also in IPECAD manuscript)
# . describe model as IPECAD1.x and IPECAD2.x


######################################## INFORMATION ########################################

# version 2.1
# see readme for details



######################################## TECHNICAL PREPARATION ########################################

cat("\014") # clear console
rm(list = ls()) # clear environment
library(dampack) # load package
# setwd("C:/users/Ron/surfdrive/PhD/Projects/IPECAD/open source model/2.1 (workshop 2023)/") # set working directory; change to your own directory
setwd("D:/surfdrive/PhD/Projects/IPECAD/open source model/2.1 (workshop 2023)/")



######################################## 1. INPUTS ########################################

######################################## 1.1. ESTIMATED ########################################

# U.S. general population life table
  # source: 
  # https://www.cdc.gov/nchs/data/nvsr/nvsr70/nvsr70-19.pdf
  # https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/NVSR/70-19/Table02.xlsx
  # https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Publications/NVSR/70-19/Table03.xlsx
m.lifetable <- as.matrix(read.csv(file="lifetable_US.csv", header=TRUE)) # import life table
m.lifetable <- m.lifetable[,c("male","female")] # select only men/women and drop age column (make sure age corresponds to row number, i.e., start with age = 1)
m.mortality_rate <- -log(1-(m.lifetable)) # convert probability to rate

# transition probability MCI to dementia in control strategy (i.e., natural disease history)
temp.p.dem18months <- 0.208 # proportion dementia as observed in IPECAD workshop control arm synthetic trial data at 18 months (source: www.ipecad.org/workshop 2023 table 3 efficacy outcomes number 'from 0.5 to 1' divided by number in trial arm = 136/654)
temp.r.dem18month <- -log(1 - temp.p.dem18months) # convert probability to rate
temp.r.dem12month <- temp.r.dem18month/1.5 # divide by 1.5 to convert 18 month interval to 12 month interval
temp.p.dem12month <- 1-exp(-temp.r.dem12month) # convert rate to probability



######################################## 1.2. MODEL INPUTS LIST ########################################

l.inputs <- list(
  v.names_state = c("mci","mil","mod","sev","dth"), # disease states: mci = mild cognitive impairment; mil = mild dementia; mod = moderate dementia; sev = severe dementia; dth = dead
  v.names_strat = c("soc","int"), # strategies: soc = standard of care strategy; int = intervention strategy
  age_start = 70, # start age (dicated by the benchmark scenario)
  age_end = 100, # end age (set at max)
  sex = "male", # sex of starting population (dependent is mortality table)
  p.mci_mil = temp.p.dem12month, # p.x_x: transition probability between states
  p.mil_mod = 0.293, 
  p.mil_sev = 0.001, 
  p.mod_mil = 0.087, 
  p.mod_sev = 0.109, 
  p.sev_mil = 0.000, 
  p.sev_mod = 0.196, 
  m.r.mortality = m.mortality_rate, # general population mortality rate
  hr.mort_mci = 1, # hr.mort_x: hazard ratio mortality by disease state (hazard ratio by dementia state compared to very mild dementia [Wimo, 2020: https://doi.org/10.3233/jad-191055] multiplied with HR or very mild compared to no dementia [Andersen, 2010: https://doi.org/10.1159/000265553])
  hr.mort_mil = 1.318 * 1.82, 
  hr.mort_mod = 2.419 * 1.82, 
  hr.mort_sev = 4.267 * 1.82, 
  rr.tx_mci_mil = 0.933, # treatment effect expressed as hazard ratio for specific transition (ratio between probability of dementia onset in control arm 0.208 and dementia onset in intervention arm 0.1942 of IPECAD workshop benchmark scenario) 0.194/0.208
  rr.tx_mil_mod = 0.933, # assumed same effect as for MCI to dementia
  p.starting_state_mci = 1, # proportion starting in disease state MCI
  u.mci = 0.73, # u.x: utility in state [https://doi.org/10.1016/j.jalz.2019.05.004]
  u.mil = 0.69, 
  u.mod = 0.53, 
  u.sev = 0.38, 
  c.mci = (1254 +  222) * 12 * (1-0    ) + (1254 + 8762) * 12 * 0, # c.x: costs in state, build up as montly costs in patient health and social care by care setting (community/residential) multiplied by 12 (annual costs) [https://doi.org/10.1007/s40120-023-00460-1] and multiplied by proportion in setting [https://doi.org/10.1007/s40120-021-00273-0]
  c.mil = (1471 +  410) * 12 * (1-0.038) + (1471 + 8762) * 12 * 0.038, 
  c.mod = (1958 +  653) * 12 * (1-0.110) + (1958 + 8762) * 12 * 0.110, 
  c.sev = (2250 + 1095) * 12 * (1-0.259) + (2250 + 8762) * 12 * 0.259, 
  c.Tx = 0, # treatment costs set to 0 according to www.ipecad.org/workshop 2023 benchmark scenario
  discount_QALY = 0.035, # discount rate set according to www.ipecad.org/workshop 2023 benchmark scenario
  discount_COST = 0.035 # # discount rate set according to www.ipecad.org/workshop 2023 benchmark scenario
)



######################################## 2. RUN MODEL ########################################

# The model is run using 2 functions: run strategy and run scenario. The second (scenario) includes a loop over all strategies by calling them one by one. This is done in the following steps (annotated in the code): 
# A: function to run a scenario
# B: prepare and initialize objects to store scenario and strategy outcomes
# C: run each strategy in a loop
# D: prepare inputs to be used in each strategy
# E: run preparations specific for the intervention strategy
# F: store newly created or updated inputs to be used for the function to run a single strategy
# G: run the strategy
# G1: prepare transition probability matrix
# G2: some checks
# G3: initialize objects to store strategy outcomes
# G4: starting state
# G5: markov multiplication by looping over cycles
# G6: multiply states with utility and cost estimates
# G7: apply discounting
# G8: apply half-cycle correction
# G9: store outcomes to be wrapped up by the 'run scenario' function
# H: store strategy results
# I: add strategy results to scenario outcomes



######################################## 2.1. RUN STRATEGY (STEP G1-G8) ########################################

# run strategy (STEP G: function for running a strategy)
f.run_strategy <- function(l.inputs) {
  with(as.list(l.inputs), {
    
    # initialize time-dependent TP matrix (STEP G1: prepare transition probability matrix)
    a.TP <- array(data = 0, dim = c(n.state, n.state, n.cycle), dimnames = list(v.names_state,v.names_state,NULL))
    
    # TP matrix state: death
    a.TP["mci","dth",] <- 1-exp(-(v.r.dth * hr.mort_mci))
    a.TP["mil","dth",] <- 1-exp(-(v.r.dth * hr.mort_mil))
    a.TP["mod","dth",] <- 1-exp(-(v.r.dth * hr.mort_mod))
    a.TP["sev","dth",] <- 1-exp(-(v.r.dth * hr.mort_sev))
    a.TP["dth","dth",] <- 1

    # TP matrix state: mci
    a.TP["mci","mci",] <- v.p.mci_mci * (1-a.TP["mci","dth",])
    a.TP["mci","mil",] <- v.p.mci_mil * (1-a.TP["mci","dth",])
    
    # TP matrix state: mild
    a.TP["mil","mil",] <- v.p.mil_mil * (1-a.TP["mil","dth",])
    a.TP["mil","mod",] <- v.p.mil_mod * (1-a.TP["mil","dth",])
    a.TP["mil","sev",] <- v.p.mil_sev * (1-a.TP["mil","dth",])
    
    # TP matrix state: moderate
    a.TP["mod","mil",] <- v.p.mod_mil * (1-a.TP["mod","dth",])
    a.TP["mod","mod",] <- v.p.mod_mod * (1-a.TP["mod","dth",])
    a.TP["mod","sev",] <- v.p.mod_sev * (1-a.TP["mod","dth",])
    
    # TP matrix state: severe
    a.TP["sev","mil",] <- v.p.sev_mil * (1-a.TP["sev","dth",])
    a.TP["sev","mod",] <- v.p.sev_mod * (1-a.TP["sev","dth",])
    a.TP["sev","sev",] <- v.p.sev_sev * (1-a.TP["sev","dth",])
    
    # check TPs sum to 1 for each cycle (STEP G2: some checks)
    for(i in v.names_state) {
      temp1 <- colSums(a.TP[i,,])
      if(!isTRUE(all.equal(current = temp1, target = rep(1,n.cycle), tolerance = 1e-10))) stop(paste("TP for",i,"do not add up to 1"))
    }
    
    # initialize output table (STEP G3: initialize objects to store strategy outcomes)
    m.out <- matrix(data = NA, nrow = n.cycle, ncol = 3, dimnames = list(NULL, c("ly","qaly","cost")))
    
    # initialize state trace
    m.trace <- matrix(data = NA, nrow = n.cycle, ncol = n.state, dimnames = list(NULL,v.names_state))
    
    # set starting state distribution (STEP G4: starting state)
    m.trace[1,] <- 0
    m.trace[1,"mci"] <- p.starting_state_mci # !!TO-DO: correct this
    m.trace[1,"mil"] <- 1-p.starting_state_mci # !!TO-DO: correct this
    
    # markov multiplication (STEP G5: markov multiplication by looping over cycles)
    for(t in 1:(n.cycle-1)) {
      m.trace[t+1,] <- m.trace[t,] %*% a.TP[,,t]
    }
    
    # primary economic outputs (STEP G6: multiply states with utility and cost estimates)
    m.out[,"ly"] <- m.trace %*% c(1,1,1,1,0)
    m.out[,"qaly"] <- m.trace %*% c(u.mci, u.mil, u.mod, u.sev, 0)
    m.out[,"cost"] <- m.trace %*% c(c.mci, c.mil, c.mod, c.sev, 0)
    
    # add additional $1000 costs for adverse events (double as this will be half-cycle corrected)
    if(strat=="int") {
      m.out[,"cost"][1] <- m.out[,"cost"][1] + 2000
    }
    
    # define vector for discounting QALYs and costs (STEP G7: apply discounting)
    v.discount_QALY <- 1 / (( 1 + (discount_QALY)) ^ (0 : (age_end-age_start-1)))
    v.discount_COST <- 1 / (( 1 + (discount_COST)) ^ (0 : (age_end-age_start-1)))
    
    # apply discounting
    m.out[,"qaly"] <- m.out[,"qaly"]*v.discount_QALY
    m.out[,"cost"] <- m.out[,"cost"]*v.discount_COST
    
    # half-cycle correction (STEP G8: apply half-cycle correction)
    for (i in 1:(n.cycle-1)) {
      m.out[i,"ly"] <- (m.out[i,"ly"]+m.out[i+1,"ly"])*0.5
      m.out[i,"qaly"] <- (m.out[i,"qaly"]+m.out[i+1,"qaly"])*0.5
      m.out[i,"cost"] <- (m.out[i,"cost"]+m.out[i+1,"cost"])*0.5
    }
    
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
f.run_scenario <- function(l.inputs) {
  with(as.list(l.inputs), { # the 'with' functions enables to call the items from the list without having to refer to the list each time (one can use 'age_start' instead l.inputs[["age_start"]])
    
    # store counters (STEP B: prepare and initialize objects to store scenario and strategy outcomes)
    n.state <- length(v.names_state) # number of states
    n.strat <- length(v.names_strat) # number of strategies
    n.cycle <- (age_end-age_start) # number of cycles
    
    # initialize output matrix (create an empty dataframe to store outcomes of a scenario)
    df.out_sum <- data.frame(
      strategy = v.names_strat,
      Cost = numeric(n.strat),
      QALY = numeric(n.strat),
      ly = numeric(n.strat),
      row.names = v.names_strat, 
      stringsAsFactors = FALSE
    )
    
    # initialize output list (create empty list to store outcomes of each strategy)
    l.out_strategy <- vector(mode = "list", length = 0)
    
    # loop over strategies (STEP C: run each strategy in a loop)
    for(strat in v.names_strat) {
      
      # convert time-independent transitions to vector of transitions (STEP D: prepare inputs to be used in each strategy)
      v.p.mci_mil <- rep(p.mci_mil, n.cycle)
      v.p.mil_mod <- rep(p.mil_mod, n.cycle)
      v.p.mil_sev <- rep(p.mil_sev, n.cycle)
      v.p.mod_mil <- rep(p.mod_mil, n.cycle)
      v.p.mod_sev <- rep(p.mod_sev, n.cycle)
      v.p.sev_mil <- rep(p.sev_mil, n.cycle)
      v.p.sev_mod <- rep(p.sev_mod, n.cycle)
      
      # probability of remaining in the same state (calculated as 1 minus transitions to other states, conditional on survival)
      v.p.mci_mci <- 1 - v.p.mci_mil
      v.p.mil_mil <- 1 - v.p.mil_mod - v.p.mil_sev
      v.p.mod_mod <- 1 - v.p.mod_mil - v.p.mod_sev
      v.p.sev_sev <- 1 - v.p.sev_mil - v.p.sev_mod
      
      # death (subset mortality table to obtain age- and sex-specific mortality)
      v.r.dth <- m.r.mortality[age_start:(age_end-1), sex]
      
      # strategy-specific inputs (STEP E: run preparations specific for the intervention strategy)
      if(strat=="int") {
        
        # waning and discontinuation
        temp.waning <- (1-0.05)^(0:(n.cycle-1)) # waning according to the www.ipecad.org/workshop 2023 benchmark scenario (set at 5% per year)
        temp.continuation <- c(1,rep(0.9,n.cycle-1)) # discontinuation according to the www.ipecad.org/workshop 2023 benchmark scenario (set at 10% after 1 year)
        temp.rr.tx <- rr.tx_mci_mil^(temp.waning * temp.continuation) # treatment effect expressed as a relative risk including a reflection of waning and discontinuation
        
        # update transition probabilities
        v.p.mci_mil <- 1-exp(-(-log(1-p.mci_mil) * temp.rr.tx)) # convert probability to rate, then multiply with treatment relative risk, then convert to probability
        v.p.mil_mod <- 1-exp(-(-log(1-p.mil_mod) * temp.rr.tx)) # idem

        # update transition probabilities of remaining in the same state
        v.p.mci_mci <- 1-v.p.mci_mil
        v.p.mil_mil <- 1-v.p.mil_mod - v.p.mil_sev
      }
      
      # list inputs for running each strategy (STEP F: store newly created or updated inputs to be used for the function to run a single strategy)
      l.inputs_strategy <- c(l.inputs, list(
          strat = strat, 
          n.state = n.state, 
          n.strat = n.strat, 
          n.cycle = n.cycle, 
          v.p.mci_mil = v.p.mci_mil, 
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
          v.r.dth = v.r.dth
      ))
      
      # run strategy (STEP G: run the strategy)
      l.strat <- f.run_strategy(l.inputs_strategy)
      
      # store strategy-specific output (STEP H store strategy results)
      l.out_strategy[[strat]] <- l.strat
      
      # store output (STEP I add strategy results to scenario outcomes)
      m.out <- l.strat[["m.out"]] # extract cycle-specific outcomes
      df.out_sum[strat,"strategy"] <- strat # store strategy name
      df.out_sum[strat,"Cost"] <- sum(m.out[,"cost"]) # calculate total costs and store them
      df.out_sum[strat,"QALY"] <- sum(m.out[,"qaly"]) # calculate total QALYs and store them
      df.out_sum[strat,"ly"] <- sum(m.out[,"ly"]) # calculate total live years and store them
      
      }
    
    # return result
    return(list(
      l.out_strategy = l.out_strategy,
      df.out_sum = df.out_sum
      )
    )
  }
)
}



######################################## 5. ANALYSIS ########################################

######################################## 5.1. PROBABILISTIC ########################################
# not applicable

######################################## 5.2. DETERMINISTIC ########################################

# run the model
out_base <- f.run_scenario(l.inputs = l.inputs)

# show output
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

# selected outcomes
round(out_base[["l.out_strategy"]][["soc" ]][["m.trace"]],4)
round(out_base[["l.out_strategy"]][["int"]][["m.trace"]],4)
out_base[["df.out_sum"]]
out_base[["df.out_sum"]][2,-1] - out_base[["df.out_sum"]][1,-1]

# icer
icer <- calculate_icers(
  cost = out_base[["df.out_sum"]][,"Cost"], 
  effect = out_base[["df.out_sum"]][,"QALY"], 
  strategies = out_base[["df.out_sum"]][,"strategy"]
)
icer

# icer plot
# plot(icer, label="all")



######################################## 5.2. DETERMINISTIC: WORKSHOP RESULTS ########################################

# results for IPECAD workshop cross-comparison
l.inputs_m <- l.inputs_f <- l.inputs
l.inputs_m[["sex"]] <- "male"
l.inputs_f[["sex"]] <- "female"
out_m <- f.run_scenario(l.inputs = l.inputs_m)
out_f <- f.run_scenario(l.inputs = l.inputs_f)

# male
m.trace_con_m <- round(out_m[["l.out_strategy"]][["soc"]][["m.trace"]],3)
m.trace_int_m <- round(out_m[["l.out_strategy"]][["int"]][["m.trace"]],3)
write.table(m.trace_con_m, "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)
write.table(m.trace_int_m, "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)
sum(out_m[["l.out_strategy"]][["soc"]][["m.out"]][1:10,"cost"])
sum(out_m[["l.out_strategy"]][["int"]][["m.out"]][1:10,"cost"])
sum(out_m[["l.out_strategy"]][["soc"]][["m.out"]][1:10,"qaly"])
sum(out_m[["l.out_strategy"]][["int"]][["m.out"]][1:10,"qaly"])

# female
m.trace_con_f <- round(out_f[["l.out_strategy"]][["soc"]][["m.trace"]],3)
m.trace_int_f <- round(out_f[["l.out_strategy"]][["int"]][["m.trace"]],3)
write.table(m.trace_con_f, "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)
write.table(m.trace_int_f, "clipboard", sep="\t", row.names=FALSE, col.names=FALSE)
sum(out_f[["l.out_strategy"]][["soc"]][["m.out"]][1:10,"cost"])
sum(out_f[["l.out_strategy"]][["int"]][["m.out"]][1:10,"cost"])
sum(out_f[["l.out_strategy"]][["soc"]][["m.out"]][1:10,"qaly"])
sum(out_f[["l.out_strategy"]][["int"]][["m.out"]][1:10,"qaly"])


