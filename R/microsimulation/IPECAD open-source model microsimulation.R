######################################## MANUAL PREPARATION ########################################

# set working directory
setwd("~/GitHub/IPECAD") # if needed, change to the directory to the folder in which the R code and the life table folder is located



######################################## TECHNICAL PREPARATION ########################################

cat("\014") # clear console
rm(list = ls()) # clear environment
gc() # garbage collection (i.e., clean up memory)



###################################### 1. INPUTS ########################################

# All model input parameters are defined and stored in a list. Some data manipulations are performed (i.e., '1.1. pre-model data manipulation') before entered in the list (i.e. '1.2. model inputs list'). 



##################################### 1.1. PRE-MODEL DATA MANIPULATION ########################################

# Prepare life table by converting probabilities to rates and put it in a convenient format before adding to the model inputs list. 

# U.S. general population life table 2019 from ssa.gov
## import life table and select only men/women and drop age column (make sure age corresponds to row number, i.e., start with age = 1)
m.lifetable_US_2019 <- as.matrix(read.csv("data/lifetable_US_2019_ssa.csv", header=TRUE))[,c("male","female")]
## convert probability to rate
m.mortality_rate_US_2019 <- -log(1-(m.lifetable_US_2019))
## weight rate for male and female
m.mortality_rate_US_2019 <- cbind(m.mortality_rate_US_2019, weighted=NA)
m.mortality_rate_US_2019[,"weighted"] <- m.mortality_rate_US_2019[,"male"] * 0.48 + m.mortality_rate_US_2019[,"female"] * 0.52
## create new mortality table with male/female
m.mortality_rate_US_2019b <- m.mortality_rate_US_2019[,c(1,2)]
## change column names to 1/2 (to match the way male/female is defined in the microsimulation)
colnames(m.mortality_rate_US_2019b) <- c(1,2)
rownames(m.mortality_rate_US_2019b) <- c(1:119)

# !!TEMPORARY: replace male/female with weighted mortality (to match IPECAD_icer Markov model replication); remove this code to use sex-specific mortality. 
m.mortality_rate_US_2019b[,1] <- m.mortality_rate_US_2019[,"weighted"]
m.mortality_rate_US_2019b[,2] <- m.mortality_rate_US_2019[,"weighted"]



######################################## 1.2. MODEL INPUTS LIST ########################################

# describe values of categorical attributes
# v.ALIVE_val <- c(0,1) # 0 = dead, 1 = alive
# v.SEX_val <- c(1,2) # 1 = male, 2 = female
# v.TX_val <- c(0,1) # 0 = intervention inactive (never provided or stopped); 1 = intervention active
# v.INSTIT_val <- c(0,1) # 0 = not living in an institution (living at home); 1 = living in an institution
# v.STATE_val <- c(1,2,3,4) # 1 = MCI, 2 = mild, 3 = moderate, 4 = severe

# input parameters
l.inputs_icer <- list(
  v.attr_names = c("STRAT","TIME","ALIVE","AGE","SEX","TX","WANING","STATE","INSTIT","QALY_PT","QALY_IC","COST_HC","COST_SC","COST_IC"),
  scenario = "ipecad icer replication", # name of the scenario
  v.strategy = c(1,2), # strategies: 1 = standard of care (soc); 2 = intervention (int)
  seed_stochastic = 12345, # seed for generating random values that drive stochastic parameters (e.g., determines baseline sex of a specific individual, or whether the event of death occurs for a certain individual at a certain time)
  seed_pa = 12345, # seed for generating random values that drive probabilistic analysis
  n.ind = 100000, # number of individuals
  n.cycle = 50, # number of cycles
  cycletime = 1, # cycle time in years
  
  # starting population
  AGE0_mean = 71, # starting population age mean
  AGE0_sd = 0, # starting population age standard deviation
  v.p.SEX0 = c(0.5,0.5), # starting population probability sex
  v.p.STATE0 = c(1,0,0,0), # starting population probability STATE
  
  # transition probabilities and relative risks
  p.mci_mil = 0.23, 
  p.mci_mod = 0, 
  p.mci_sev = 0, 
  p.mil_mci = 0.03, 
  p.mil_mod = 0.35, 
  p.mil_sev = 0.04, 
  p.mod_mci = 0, 
  p.mod_mil = 0.03, 
  p.mod_sev = 0.42, 
  p.sev_mci = 0, 
  p.sev_mil = 0, 
  p.sev_mod = 0.02, 
  p.mci_i = 0.024, 
  p.mil_i = 0.038, 
  p.mod_i = 0.110, 
  p.sev_i = 0.259, 
  m.lifetable = m.mortality_rate_US_2019b,
  v.hr_mort_state = c(1.82, 2.92, 3.85, 9.52), 
  
  # treatment effect
  rr.Tx_mci_mil = 0.69, 
  rr.Tx_mci_mod = 1, 
  rr.Tx_mci_sev = 1, 
  rr.Tx_mil_mod = 0.69, 
  rr.Tx_mil_sev = 0.69, 
  rr.Tx_mci_mil_dis = 1, 
  rr.Tx_mci_mod_dis = 1, 
  rr.Tx_mci_sev_dis = 1, 
  rr.Tx_mil_mod_dis = 1, 
  rr.Tx_mil_sev_dis = 1, 
  p.tx_discontinuation1 = 0.2,  
  p.tx_discontinuation2 = 0.5,  
  tx_discontinuation2_begin = 2, # !!TO-DO: this feature needs to be programmed into the model function
  tx_duration = 10, 
  tx_waning = 0.8, 
  tx_waning_dis = 0.2, 
  
  # utilities
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
  
  # costs
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
  c.Tx = 0, # 26500 + (52/2)*78.35, # drug annual wholesale acquisition cost + treatment administration frequency * administration cost
  c.Tx_start = 261.10*4 + 261.10*3*0.215, # mri cost * 3-month monitoring in year 1 + mri cost * 3 times * proportion aria
  
  # discount rate and willingness to pay
  discount_EFFECT = 0, # !!TO-DO: this feature needs to be programmed into the model function
  discount_QALY = 0, # !!TO-DO: this feature needs to be programmed into the model function
  discount_COST = 0, # !!TO-DO: this feature needs to be programmed into the model function
  wtp = 10000
)

# add names to objects
names(l.inputs_icer[["v.p.SEX0"]]) <- c("1","2")
names(l.inputs_icer[["v.hr_mort_state"]]) <- c("1","2","3","4")
names(l.inputs_icer[["v.p.STATE0"]]) <- c("mci","mil","mod","sev")



######################################## 2. RUN MODEL ########################################

# Main section for running the microsimulation model.



###################################### 2.1.1. FUNCTIONS ########################################

# Define generic and model-specific functions used throughout the simulation.

# probability to rate
f.pr <- function(p, t=1) {
  if ((sum(p > 1, na.rm=TRUE) > 0) | (sum(p < 0, na.rm=TRUE) > 0)) stop("probability not between 0 and 1")
  r = -(1/t)*log(1 - p)
  return(r)
}

# rate to probability
f.rp <- function(r, t=1) {
  if ((sum(r < 0, na.rm=TRUE) > 0)) stop("rate not greater than or equal to 0")
  p <- 1 - exp(-r * t)
  return(p)
}

# function: discounting
f.discount <- function(x, discount_rate, n.cycle) {
  as.matrix(x) / (1 + discount_rate)^(0:(n.cycle - 1)) # x can be scalar, vector or matrix
}

# function: draw random values for each attribute for each cycle for each individual
f.random <- function(n.ind, n.cycle, n.attr, seed_stochastic, v.attr_names) {
  a.out <- array(
    data = NA,
    dim = c(n.cycle, n.attr, n.ind),
    dimnames = list(NULL, v.attr_names, NULL)
  )
  set.seed(seed_stochastic)
  a.out[,,] <- runif(n = length(a.out))
  return(a.out)
}



######################################## 2.1.2. RUN SCENARIO ########################################

# Main function to run a scenario for all strategies, cycles, and individuals, and collect results.

# function to run a scenario
f.run_scenario <- function(l.inputs, detailed=FALSE) {
  with(as.list(l.inputs), {
    
    # set start time
    ptm <- proc.time()
    stime <- Sys.time()
    
    # random values at each cycle for each attribute for each individual and store in array called 'a.rnd'
    n.attr <- length(v.attr_names)
    a.rnd <- f.random(n.ind=n.ind, n.cycle=n.cycle, n.attr=n.attr, seed_stochastic=seed_stochastic, v.attr_names=v.attr_names)
    
    # initialize: data frame for outcomes of each strategy (STEP A1: initialize objects to store scenario and strategy outcomes)
    v.names_strat <- v.strategy
    n.strat <- length(v.names_strat) # number of strategies
    df.out <- data.frame(
      strategy = v.names_strat,
      QALY = numeric(n.strat),
      COST = numeric(n.strat),
      LY = numeric(n.strat),
      NHB = numeric(n.strat),
      row.names = v.names_strat,
      stringsAsFactors = FALSE
    )
    
    # initialize: list for outcomes of each scenario
    l.out_scenario <- vector(mode = "list", length = 0)
    
    # loop over strategies (STEP A2: run each strategy in a loop)
    for(strategy in v.strategy) {
      
      # run strategy (STEP B: run the strategy)
      # !!TO-DO: from this point onwards the steps should run in a separate function 'f.run_strategy' to enable the functionality of 'dampack' package, see IPECAD markov model
      # l.out_strategy <- f.run_strategy(l.inputs = l.inputs, strat = strat)
      
      # initialize: array for outcomes of each individual
      a.out <- a.rnd
      a.out[,,] <- NA
      
      # within-model data manipulation
      ## transition probability matrix: soc
      m.TP <- matrix(
        data = c(
          NA       , p.mci_mil, p.mci_mod, p.mci_sev, 
          p.mil_mci, NA       , p.mil_mod, p.mil_sev, 
          p.mod_mci, p.mod_mil, NA       , p.mod_sev, 
          p.sev_mci, p.sev_mil, p.sev_mod, NA       
        ),
        byrow = TRUE, 
        nrow = 4, 
        ncol = 4,
        dimnames = list(c("mci","mil","mod","sev"),c("mci","mil","mod","sev")) # transitions reflect from row to column
      )
      ## matrix: relative risk treatment effect
      m.RR_tx <- matrix(
        data = c(
          1            , rr.Tx_mci_mil, rr.Tx_mci_mod, rr.Tx_mci_sev, 
          1            , 1            , rr.Tx_mil_mod, rr.Tx_mil_sev, 
          1            , 1            , 1            , 1            , 
          1            , 1            , 1            , 1
        ),
        byrow = TRUE, 
        nrow = 4, 
        ncol = 4,
        dimnames = list(c("mci","mil","mod","sev"),c("mci","mil","mod","sev")) # transitions reflect from row to column
      )
      ## matrix: relative risk treatment effect after discontinuation
      m.RR_tx_dis <- matrix(
        data = c(
          1            , rr.Tx_mci_mil_dis, rr.Tx_mci_mod_dis, rr.Tx_mci_sev_dis, 
          1            , 1                , rr.Tx_mil_mod_dis, rr.Tx_mil_sev_dis, 
          1            , 1                , 1                , 1                , 
          1            , 1                , 1                , 1
        ),
        byrow = TRUE, 
        nrow = 4, 
        ncol = 4,
        dimnames = list(c("mci","mil","mod","sev"),c("mci","mil","mod","sev")) # transitions reflect from row to column
      )
      
      # record strategy (all cycles simultaneously)
      a.out[,"STRAT",] <- strategy
      
      # record time (all cycles simultaneously)
      a.out[,"TIME",] <- rep(x=1:n.cycle, times=n.ind)
      
      # set the baseline values for all individuals in the first cycle
      a.out[1,"ALIVE",] <- 1
      a.out[1,"AGE",] <- AGE0_mean
      a.out[1,"SEX",] <- .bincode(x=a.rnd[1,"SEX",], breaks=c(0, cumsum(v.p.SEX0)), include.lowest=TRUE)
      a.out[1,"TX",] <- ifelse(test=strategy==2, yes=1, no=0)
        #a.out[1,"WANING",] <- ifelse(test=strategy==2, yes=0.10, no=1) # !!TO-DO: this is not yet implemented in this microsimulation
      a.out[1,"STATE",] <- sample(c(1,2,3,4), size = n.ind, replace = TRUE, prob = v.p.STATE0)
      a.out[1,"INSTIT",] <- 0
      
      # set baseline values for QALY_PT, QALY_IC, COST_HC, COST_SC, COST_IC
      a.out[1,"QALY_PT",] <- 
        (a.out[1,"STATE",]==1) * ((a.out[1,"INSTIT",]==0) * u.mci_pt + (a.out[1,"INSTIT",]==1) * u.mci_pt_i) +
        (a.out[1,"STATE",]==2) * ((a.out[1,"INSTIT",]==0) * u.mil_pt + (a.out[1,"INSTIT",]==1) * u.mil_pt_i) +
        (a.out[1,"STATE",]==3) * ((a.out[1,"INSTIT",]==0) * u.mod_pt + (a.out[1,"INSTIT",]==1) * u.mod_pt_i) +
        (a.out[1,"STATE",]==4) * ((a.out[1,"INSTIT",]==0) * u.sev_pt + (a.out[1,"INSTIT",]==1) * u.sev_pt_i) +
        u.Tx_start
      a.out[1,"QALY_IC",] <-
        (a.out[1,"STATE",]==1) * ((a.out[1,"INSTIT",]==0) * u.mci_ic + (a.out[1,"INSTIT",]==1) * u.mci_ic_i) +
        (a.out[1,"STATE",]==2) * ((a.out[1,"INSTIT",]==0) * u.mil_ic + (a.out[1,"INSTIT",]==1) * u.mil_ic_i) +
        (a.out[1,"STATE",]==3) * ((a.out[1,"INSTIT",]==0) * u.mod_ic + (a.out[1,"INSTIT",]==1) * u.mod_ic_i) +
        (a.out[1,"STATE",]==4) * ((a.out[1,"INSTIT",]==0) * u.sev_ic + (a.out[1,"INSTIT",]==1) * u.sev_ic_i)
      a.out[1,"COST_HC",] <-
        (a.out[1,"STATE",]==1) * ((a.out[1,"INSTIT",]==0) * c.mci_hc + (a.out[1,"INSTIT",]==1) * c.mci_hc_i) +
        (a.out[1,"STATE",]==2) * ((a.out[1,"INSTIT",]==0) * c.mil_hc + (a.out[1,"INSTIT",]==1) * c.mil_hc_i) +
        (a.out[1,"STATE",]==3) * ((a.out[1,"INSTIT",]==0) * c.mod_hc + (a.out[1,"INSTIT",]==1) * c.mod_hc_i) +
        (a.out[1,"STATE",]==4) * ((a.out[1,"INSTIT",]==0) * c.sev_hc + (a.out[1,"INSTIT",]==1) * c.sev_hc_i) + 
        c.Tx_start
      a.out[1,"COST_SC",] <-
        (a.out[1,"STATE",]==1) * ((a.out[1,"INSTIT",]==0) * c.mci_sc + (a.out[1,"INSTIT",]==1) * c.mci_sc_i) +
        (a.out[1,"STATE",]==2) * ((a.out[1,"INSTIT",]==0) * c.mil_sc + (a.out[1,"INSTIT",]==1) * c.mil_sc_i) +
        (a.out[1,"STATE",]==3) * ((a.out[1,"INSTIT",]==0) * c.mod_sc + (a.out[1,"INSTIT",]==1) * c.mod_sc_i) +
        (a.out[1,"STATE",]==4) * ((a.out[1,"INSTIT",]==0) * c.sev_sc + (a.out[1,"INSTIT",]==1) * c.sev_sc_i)
      a.out[1,"COST_IC",] <-
        (a.out[1,"STATE",]==1) * ((a.out[1,"INSTIT",]==0) * c.mci_ic + (a.out[1,"INSTIT",]==1) * c.mci_ic_i) +
        (a.out[1,"STATE",]==2) * ((a.out[1,"INSTIT",]==0) * c.mil_ic + (a.out[1,"INSTIT",]==1) * c.mil_ic_i) +
        (a.out[1,"STATE",]==3) * ((a.out[1,"INSTIT",]==0) * c.mod_ic + (a.out[1,"INSTIT",]==1) * c.mod_ic_i) +
        (a.out[1,"STATE",]==4) * ((a.out[1,"INSTIT",]==0) * c.sev_ic + (a.out[1,"INSTIT",]==1) * c.sev_ic_i)
      
      # progress bar
      pb = txtProgressBar(min = 0, max = n.cycle * length(v.strategy), initial = 0) # set progress bar
      setTxtProgressBar(pb, 1) # update progress bar
      
      # loop over cycles
      for (i in 2:n.cycle) {
        
        # select individuals: alive at previous cycle
        v.alive.lag <- a.out[i-1,"ALIVE",]== 1
        
        # update: ALIVE (if dead at previous cycle)
        a.out[i,"ALIVE",a.out[i-1,"ALIVE",]==0] <- 0
        
        # update: ALIVE
        m.lookupcoordinates_lifetable <- matrix(data = as.character(c(a.out[i-1,"AGE",v.alive.lag], a.out[i-1,"SEX",v.alive.lag])), ncol = 2) # generate life table coordinates for looking up age- and sex-specific mortality; (see https://adv-r.hadley.nz/subsetting.html paragraph 4.2.3 subsetting > selecting multiple elements > subsetting)
        v.lookupcoordinates_RR <- as.character(a.out[i-1,"STATE",v.alive.lag]) # determine relative mortality risk related to syndrome and severity
        v.p_dth <- 1 - exp(-m.lifetable[m.lookupcoordinates_lifetable] * v.hr_mort_state[v.lookupcoordinates_RR]) # probability of death (input = lifetable rate adjusted for HR dementia severity)
        a.out[i,"ALIVE",v.alive.lag] <- as.numeric(!v.p_dth > a.rnd[i,"ALIVE",v.alive.lag]) # compare probability to random value
        
        # select individuals: alive at current cycle
        v.alive <- a.out[i,"ALIVE",]==1
        
        # update: AGE
        a.out[i,"AGE",v.alive] <- a.out[i-1,"AGE",v.alive] + cycletime
        
        # update: SEX
        a.out[i,"SEX",v.alive] <- a.out[i-1,"SEX",v.alive]
        
        # update: TX
        if(strategy==1) a.out[i,"TX",v.alive] <- 0
        if(strategy==2) {
          a.out[i,"TX",v.alive] <- 
            a.out[i-1,"TX",v.alive]==1 & 
            ( (a.out[i-1,"TIME",v.alive]==1 & a.rnd[i,"TX",v.alive]>p.tx_discontinuation1)*1 + (a.out[i-1,"TIME",v.alive]>1 & a.rnd[i,"TX",v.alive]>p.tx_discontinuation2)*1 ) & 
            a.out[i-1,"INSTIT",v.alive]==0 & 
            a.out[i-1,"TIME",v.alive]<tx_duration & 
            (a.out[i-1,"STATE",v.alive]==1 | a.out[i-1,"STATE",v.alive]==2)
        }
        
        # update: STATE
        # alternative could be to individualize all 16 transitions and track them in separate columns as an attribute of the individual, or program this in an individual loop (rather than to vectorize this)
        ## transition probability matrix
        m.TP_soc    <- m.TP # soc
        m.TP_tx     <- f.rp(f.pr(m.TP_soc) * m.RR_tx^((1-tx_waning)^(i-2))) # int: treatment
        m.TP_tx_dis <- f.rp(f.pr(m.TP_soc) * m.RR_tx_dis^((1-tx_waning_dis)^(i-2))) # int: treatment discontinued
        ## remain in same state
        diag(m.TP_soc) <- 1 - rowSums(m.TP_soc, na.rm=TRUE)
        diag(m.TP_tx) <- 1 - rowSums(m.TP_tx, na.rm=TRUE)
        diag(m.TP_tx_dis) <- 1 - rowSums(m.TP_tx_dis, na.rm=TRUE)
        ## cumulative probabilities from each state
        m.TP_soc_breaks    <- cbind(0, t(apply(X=m.TP_soc   , MARGIN=1, FUN=cumsum))) # cumulative probabilities for each 'from' state
        m.TP_tx_breaks     <- cbind(0, t(apply(X=m.TP_tx    , MARGIN=1, FUN=cumsum))) # cumulative probabilities for each 'from' state
        m.TP_tx_dis_breaks <- cbind(0, t(apply(X=m.TP_tx_dis, MARGIN=1, FUN=cumsum))) # cumulative probabilities for each 'from' state
        ## update state
        if(strategy==1) {
          a.out[i,"STATE",v.alive & a.out[i-1,"STATE",]==1] <- .bincode(x=a.rnd[i,"STATE",v.alive & a.out[i-1,"STATE",]==1], breaks=m.TP_soc_breaks[1,], include.lowest=TRUE)
          a.out[i,"STATE",v.alive & a.out[i-1,"STATE",]==2] <- .bincode(x=a.rnd[i,"STATE",v.alive & a.out[i-1,"STATE",]==2], breaks=m.TP_soc_breaks[2,], include.lowest=TRUE)
          a.out[i,"STATE",v.alive & a.out[i-1,"STATE",]==3] <- .bincode(x=a.rnd[i,"STATE",v.alive & a.out[i-1,"STATE",]==3], breaks=m.TP_soc_breaks[3,], include.lowest=TRUE)
          a.out[i,"STATE",v.alive & a.out[i-1,"STATE",]==4] <- .bincode(x=a.rnd[i,"STATE",v.alive & a.out[i-1,"STATE",]==4], breaks=m.TP_soc_breaks[4,], include.lowest=TRUE)
          # # alternative 1
          # v.rnd <- a.rnd[i,"STATE",v.alive] # temporary vector
          # v.state.lag <- a.out[i-1,"STATE",v.alive] # temporary vector
          # v.out_state <- rep(NA, length(v.rnd)) # temporary vector
          # for(x in 1:nrow(m.TP_soc) ) v.out_state[v.state.lag==x] <- .bincode(x=v.rnd[v.state.lag==x], breaks=m.TP_soc_breaks[x,], include.lowest=TRUE)
          # a.out[i,"STATE",v.alive] <- v.out_state
          # # alternative 2 (does not work)
          # a.out[i,"STATE",v.alive] <- 
          #   (a.out[i-1,"STATE",v.alive]==1) * .bincode(x=a.rnd[i,"STATE",v.alive & a.out[i-1,"STATE",]==1], breaks=m.TP_soc_breaks[1,], include.lowest=TRUE) + 
          #   (a.out[i-1,"STATE",v.alive]==2) * .bincode(x=a.rnd[i,"STATE",v.alive & a.out[i-1,"STATE",]==2], breaks=m.TP_soc_breaks[2,], include.lowest=TRUE) + 
          #   (a.out[i-1,"STATE",v.alive]==3) * .bincode(x=a.rnd[i,"STATE",v.alive & a.out[i-1,"STATE",]==3], breaks=m.TP_soc_breaks[3,], include.lowest=TRUE) + 
          #   (a.out[i-1,"STATE",v.alive]==4) * .bincode(x=a.rnd[i,"STATE",v.alive & a.out[i-1,"STATE",]==4], breaks=m.TP_soc_breaks[4,], include.lowest=TRUE)
        } else if(strategy==2) {
          a.out[i,"STATE",v.alive & a.out[i-1,"TX",]==1 & a.out[i-1,"STATE",]==1] <- .bincode(x=a.rnd[i,"STATE",v.alive & a.out[i-1,"TX",]==1 & a.out[i-1,"STATE",]==1], breaks=m.TP_tx_breaks[1,], include.lowest=TRUE)
          a.out[i,"STATE",v.alive & a.out[i-1,"TX",]==1 & a.out[i-1,"STATE",]==2] <- .bincode(x=a.rnd[i,"STATE",v.alive & a.out[i-1,"TX",]==1 & a.out[i-1,"STATE",]==2], breaks=m.TP_tx_breaks[2,], include.lowest=TRUE)
          a.out[i,"STATE",v.alive & a.out[i-1,"TX",]==1 & a.out[i-1,"STATE",]==3] <- .bincode(x=a.rnd[i,"STATE",v.alive & a.out[i-1,"TX",]==1 & a.out[i-1,"STATE",]==3], breaks=m.TP_tx_breaks[3,], include.lowest=TRUE)
          a.out[i,"STATE",v.alive & a.out[i-1,"TX",]==1 & a.out[i-1,"STATE",]==4] <- .bincode(x=a.rnd[i,"STATE",v.alive & a.out[i-1,"TX",]==1 & a.out[i-1,"STATE",]==4], breaks=m.TP_tx_breaks[4,], include.lowest=TRUE)
          a.out[i,"STATE",v.alive & a.out[i-1,"TX",]==0 & a.out[i-1,"STATE",]==1] <- .bincode(x=a.rnd[i,"STATE",v.alive & a.out[i-1,"TX",]==0 & a.out[i-1,"STATE",]==1], breaks=m.TP_tx_dis_breaks[1,], include.lowest=TRUE)
          a.out[i,"STATE",v.alive & a.out[i-1,"TX",]==0 & a.out[i-1,"STATE",]==2] <- .bincode(x=a.rnd[i,"STATE",v.alive & a.out[i-1,"TX",]==0 & a.out[i-1,"STATE",]==2], breaks=m.TP_tx_dis_breaks[2,], include.lowest=TRUE)
          a.out[i,"STATE",v.alive & a.out[i-1,"TX",]==0 & a.out[i-1,"STATE",]==3] <- .bincode(x=a.rnd[i,"STATE",v.alive & a.out[i-1,"TX",]==0 & a.out[i-1,"STATE",]==3], breaks=m.TP_tx_dis_breaks[3,], include.lowest=TRUE)
          a.out[i,"STATE",v.alive & a.out[i-1,"TX",]==0 & a.out[i-1,"STATE",]==4] <- .bincode(x=a.rnd[i,"STATE",v.alive & a.out[i-1,"TX",]==0 & a.out[i-1,"STATE",]==4], breaks=m.TP_tx_dis_breaks[4,], include.lowest=TRUE)
        }
        
        # update: INSTIT
        a.out[i,"INSTIT",v.alive & a.out[i-1,"INSTIT",]==1] <- 1
        a.out[i,"INSTIT",v.alive & a.out[i-1,"INSTIT",]==0 & a.out[i-1,"STATE",]==1] <- ifelse(test = a.rnd[i,"INSTIT",v.alive & a.out[i-1,"INSTIT",]==0 & a.out[i-1,"STATE",]==1] < p.mci_i, yes=1, no=0)
        a.out[i,"INSTIT",v.alive & a.out[i-1,"INSTIT",]==0 & a.out[i-1,"STATE",]==2] <- ifelse(test = a.rnd[i,"INSTIT",v.alive & a.out[i-1,"INSTIT",]==0 & a.out[i-1,"STATE",]==2] < p.mil_i, yes=1, no=0)
        a.out[i,"INSTIT",v.alive & a.out[i-1,"INSTIT",]==0 & a.out[i-1,"STATE",]==3] <- ifelse(test = a.rnd[i,"INSTIT",v.alive & a.out[i-1,"INSTIT",]==0 & a.out[i-1,"STATE",]==3] < p.mod_i, yes=1, no=0)
        a.out[i,"INSTIT",v.alive & a.out[i-1,"INSTIT",]==0 & a.out[i-1,"STATE",]==4] <- ifelse(test = a.rnd[i,"INSTIT",v.alive & a.out[i-1,"INSTIT",]==0 & a.out[i-1,"STATE",]==4] < p.sev_i, yes=1, no=0)
        
        # update: QALY_PT
        a.out[i,"QALY_PT",v.alive] <- 
          (a.out[i,"STATE",v.alive]==1) * ((a.out[i,"INSTIT",v.alive]==0) * u.mci_pt + (a.out[i,"INSTIT",v.alive]==1) * u.mci_pt_i) +
          (a.out[i,"STATE",v.alive]==2) * ((a.out[i,"INSTIT",v.alive]==0) * u.mil_pt + (a.out[i,"INSTIT",v.alive]==1) * u.mil_pt_i) +
          (a.out[i,"STATE",v.alive]==3) * ((a.out[i,"INSTIT",v.alive]==0) * u.mod_pt + (a.out[i,"INSTIT",v.alive]==1) * u.mod_pt_i) +
          (a.out[i,"STATE",v.alive]==4) * ((a.out[i,"INSTIT",v.alive]==0) * u.sev_pt + (a.out[i,"INSTIT",v.alive]==1) * u.sev_pt_i)
        
        # update: QALY_IC
        a.out[i,"QALY_IC",v.alive] <- 
          (a.out[i,"STATE",v.alive]==1) * ((a.out[i,"INSTIT",v.alive]==0) * u.mci_ic + (a.out[i,"INSTIT",v.alive]==1) * u.mci_ic_i) +
          (a.out[i,"STATE",v.alive]==2) * ((a.out[i,"INSTIT",v.alive]==0) * u.mil_ic + (a.out[i,"INSTIT",v.alive]==1) * u.mil_ic_i) +
          (a.out[i,"STATE",v.alive]==3) * ((a.out[i,"INSTIT",v.alive]==0) * u.mod_ic + (a.out[i,"INSTIT",v.alive]==1) * u.mod_ic_i) +
          (a.out[i,"STATE",v.alive]==4) * ((a.out[i,"INSTIT",v.alive]==0) * u.sev_ic + (a.out[i,"INSTIT",v.alive]==1) * u.sev_ic_i)
        
        # update: COST_HC
        a.out[i,"COST_HC",v.alive] <- 
          (a.out[i,"STATE",v.alive]==1) * ((a.out[i,"INSTIT",v.alive]==0) * c.mci_hc + (a.out[i,"INSTIT",v.alive]==1) * c.mci_hc_i) +
          (a.out[i,"STATE",v.alive]==2) * ((a.out[i,"INSTIT",v.alive]==0) * c.mil_hc + (a.out[i,"INSTIT",v.alive]==1) * c.mil_hc_i) +
          (a.out[i,"STATE",v.alive]==3) * ((a.out[i,"INSTIT",v.alive]==0) * c.mod_hc + (a.out[i,"INSTIT",v.alive]==1) * c.mod_hc_i) +
          (a.out[i,"STATE",v.alive]==4) * ((a.out[i,"INSTIT",v.alive]==0) * c.sev_hc + (a.out[i,"INSTIT",v.alive]==1) * c.sev_hc_i)
        
        # update: COST_SC
        a.out[i,"COST_SC",v.alive] <- 
          (a.out[i,"STATE",v.alive]==1) * ((a.out[i,"INSTIT",v.alive]==0) * c.mci_sc + (a.out[i,"INSTIT",v.alive]==1) * c.mci_sc_i) +
          (a.out[i,"STATE",v.alive]==2) * ((a.out[i,"INSTIT",v.alive]==0) * c.mil_sc + (a.out[i,"INSTIT",v.alive]==1) * c.mil_sc_i) +
          (a.out[i,"STATE",v.alive]==3) * ((a.out[i,"INSTIT",v.alive]==0) * c.mod_sc + (a.out[i,"INSTIT",v.alive]==1) * c.mod_sc_i) +
          (a.out[i,"STATE",v.alive]==4) * ((a.out[i,"INSTIT",v.alive]==0) * c.sev_sc + (a.out[i,"INSTIT",v.alive]==1) * c.sev_sc_i)
        
        # update: COST_IC
        a.out[i,"COST_IC",v.alive] <- 
          (a.out[i,"STATE",v.alive]==1) * ((a.out[i,"INSTIT",v.alive]==0) * c.mci_ic + (a.out[i,"INSTIT",v.alive]==1) * c.mci_ic_i) +
          (a.out[i,"STATE",v.alive]==2) * ((a.out[i,"INSTIT",v.alive]==0) * c.mil_ic + (a.out[i,"INSTIT",v.alive]==1) * c.mil_ic_i) +
          (a.out[i,"STATE",v.alive]==3) * ((a.out[i,"INSTIT",v.alive]==0) * c.mod_ic + (a.out[i,"INSTIT",v.alive]==1) * c.mod_ic_i) +
          (a.out[i,"STATE",v.alive]==4) * ((a.out[i,"INSTIT",v.alive]==0) * c.sev_ic + (a.out[i,"INSTIT",v.alive]==1) * c.sev_ic_i)
        
        # update progress bar
        setTxtProgressBar(pb, i)
      }
      
      # initialize: list for strategy results
      l.out_strategy <- vector(mode = "list", length = 0)
      
      # store results strategy
      l.out_strategy[["a.out"]] <- a.out # individual patient level data
      ## trace
      m.trace <- matrix(data=NA, nrow=n.cycle, ncol=20, dimnames=list(NULL,c("mcion_c","mciof_c","milon_c","milof_c","mod_c","sev_c","mci_i","mil_i","mod_i","sev_i","dth","qaly_pt","qaly_ic","cost_hc","cost_sc","cost_ic","alv","instit","qaly_total","cost_total")))
      m.trace[,"mcion_c"] <- as.matrix(apply(X = a.out[,"STATE",]==1 & a.out[,"TX",]==1 & a.out[,"INSTIT",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n.ind)
      m.trace[,"mciof_c"] <- as.matrix(apply(X = a.out[,"STATE",]==1 & a.out[,"TX",]==0 & a.out[,"INSTIT",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n.ind)
      m.trace[,"milon_c"] <- as.matrix(apply(X = a.out[,"STATE",]==2 & a.out[,"TX",]==1 & a.out[,"INSTIT",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n.ind)
      m.trace[,"milof_c"] <- as.matrix(apply(X = a.out[,"STATE",]==2 & a.out[,"TX",]==0 & a.out[,"INSTIT",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n.ind)
      m.trace[,"mod_c"] <- as.matrix(apply(X = a.out[,"STATE",]==3 & a.out[,"INSTIT",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n.ind)
      m.trace[,"sev_c"] <- as.matrix(apply(X = a.out[,"STATE",]==4 & a.out[,"INSTIT",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n.ind)
      m.trace[,"mci_i"] <- as.matrix(apply(X = a.out[,"STATE",]==1 & a.out[,"INSTIT",]==1, MARGIN = 1, FUN = sum, na.rm = TRUE)/n.ind)
      m.trace[,"mil_i"] <- as.matrix(apply(X = a.out[,"STATE",]==2 & a.out[,"INSTIT",]==1, MARGIN = 1, FUN = sum, na.rm = TRUE)/n.ind)
      m.trace[,"mod_i"] <- as.matrix(apply(X = a.out[,"STATE",]==3 & a.out[,"INSTIT",]==1, MARGIN = 1, FUN = sum, na.rm = TRUE)/n.ind)
      m.trace[,"sev_i"] <- as.matrix(apply(X = a.out[,"STATE",]==4 & a.out[,"INSTIT",]==1, MARGIN = 1, FUN = sum, na.rm = TRUE)/n.ind)
      m.trace[,"dth"] <- as.matrix(apply(X = a.out[,"ALIVE",]==0, MARGIN = 1, FUN = sum, na.rm = TRUE)/n.ind)
      m.trace[,"qaly_pt"] <- as.matrix(apply(X = a.out[,"QALY_PT",], MARGIN = 1, FUN = sum, na.rm = TRUE)/n.ind)
      m.trace[,"qaly_ic"] <- as.matrix(apply(X = a.out[,"QALY_IC",], MARGIN = 1, FUN = sum, na.rm = TRUE)/n.ind)
      m.trace[,"cost_hc"] <- as.matrix(apply(X = a.out[,"COST_HC",], MARGIN = 1, FUN = sum, na.rm = TRUE)/n.ind)
      m.trace[,"cost_sc"] <- as.matrix(apply(X = a.out[,"COST_SC",], MARGIN = 1, FUN = sum, na.rm = TRUE)/n.ind)
      m.trace[,"cost_ic"] <- as.matrix(apply(X = a.out[,"COST_IC",], MARGIN = 1, FUN = sum, na.rm = TRUE)/n.ind)
      m.trace[,"alv"] <- as.matrix(apply(X = a.out[,"ALIVE",]==1, MARGIN = 1, FUN = sum, na.rm = TRUE)/n.ind)
      m.trace[,"instit"] <- rowSums(m.trace[,c("mci_i","mil_i","mod_i","sev_i")])
      m.trace[,"qaly_total"] <- rowSums(m.trace[,c("qaly_pt","qaly_ic")])
      m.trace[,"cost_total"] <- rowSums(m.trace[,c("cost_hc","cost_sc","cost_ic")])
      l.out_strategy[["m.trace"]] <- m.trace
      
      # store strategy-specific output (STEP A3 store strategy results)
      l.out_scenario[[strategy]] <- l.out_strategy
      
      # store output (STEP A4 add strategy results to scenario outcomes)
      df.out[strategy,"strategy"] <- strategy # store strategy name
      df.out[strategy,"QALY"] <- sum(m.trace[,"qaly_total"])
      df.out[strategy,"COST"] <- sum(m.trace[,"cost_total"])
      df.out[strategy,"LY"]   <- sum(m.trace[,"alv"])
      df.out[strategy,"NHB"]  <- df.out[strategy,"QALY"] - (df.out[strategy,"COST"] / wtp) # calculate total NHB and store them
      
    }
    
    # run time
    print(proc.time() - ptm)
    print(Sys.time() - stime)
    
    # return result
    if(detailed) return(list(df.out = df.out, l.out_scenario = l.out_scenario))
    if(!detailed) return(df.out)
    
  })
}


######################################## 2.2. RUN STRATEGY ########################################
# Placeholder for a function to run a single strategy (not yet implemented).


######################################## 3. MODEL CALIBRATION ########################################
# Section for model calibration routines (currently not applicable).


######################################## 4. VALIDATION ########################################
# Section for model validation routines (currently not applicable).


######################################## 5. ANALYSIS ########################################
# Run the base case simulation and display or save results.

# run model: base case
l.out_base <- f.run_scenario(l.inputs=l.inputs_icer, detailed=TRUE)



######################################## 6. RESULT ########################################

# base case results
l.out_base[["df.out"]]
str(l.out_base) # print structure of all data stored in the outputs list

#Individual level data for both strategies
l.out_base[["l.out_scenario"]][[1]][["a.out"]] # individual level data for strategy 1
l.out_base[["l.out_scenario"]][[2]][["a.out"]] # individual level data for strategy 2

# state trace
l.out_base[["l.out_scenario"]][[1]][["m.trace"]]
l.out_base[["l.out_scenario"]][[2]][["m.trace"]]

# internal validity check for health states
rowSums(l.out_base[["l.out_scenario"]][[1]][["m.trace"]][,c("mcion_c","mciof_c","milon_c","milof_c","mod_c","sev_c","mci_i","mil_i","mod_i","sev_i","dth")]) # check if add up to 1
rowSums(l.out_base[["l.out_scenario"]][[2]][["m.trace"]][,c("mcion_c","mciof_c","milon_c","milof_c","mod_c","sev_c","mci_i","mil_i","mod_i","sev_i","dth")]) # check if add up to 1

# aggregated outcomes for all individuals
colSums(l.out_base[["l.out_scenario"]][[1]][["m.trace"]])
colSums(l.out_base[["l.out_scenario"]][[2]][["m.trace"]])
