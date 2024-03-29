# load packages
library(shiny)

######################################## UI ########################################

ui <- fluidPage(
  
  # App title
  titlePanel("IPECAD open-source model"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    
    
    
    ######################################## * SIDEBAR PANEL ########################################
    
    sidebarPanel(
    
      h3("Basic inputs"), 
      
      helpText(
        ("Abbreviations:"), br(), 
        ("HR = hazard ratio"), br(),
        ("MCI = mild cognitive impairment"), br(),
        ("RR = relative risk"),
        ("TP = transition probability"), br(), 
      ), 
      
      sliderInput(inputId = "age_start", label = "Age of starting population", min = 50, max = 99, value = 71), 
      radioButtons(inputId = "sex", label = "Sex of starting population", choices = c("female","male","weighted"), selected = "weighted"), 
      sliderInput(inputId = "p.starting_state_mci", label="Proportion starting population in MCI", min = 0, max = 1, value = 0.55), 
      sliderInput(inputId = "n.cycle", label="Number of cycles (years) to run", min = 1, max = 50, value = 29), 
      sliderInput(inputId = "p.mci_mil", label = "TP MCI to mild dementia", min = 0, max = 1, value = 0.23), 
      sliderInput(inputId = "p.mil_mod", label = "TP mild to moderate dementia", min = 0, max = 1, value = 0.35), 
      sliderInput(inputId = "rr.tx_mci_mil", label = "RR treatment effect MCI to mild dementia", min = 0, max = 1, value = 0.69), 
      sliderInput(inputId = "rr.tx_mil_mod", label = "RR treatment effect mild dementia to moderate dementia", min = 0, max = 1, value = 0.69), 
      sliderInput(inputId = "rr.tx_mil_sev", label = "RR treatment effect mild dementia to severe dementia", min = 0, max = 1, value = 0.69), 
      numericInput(inputId = "p.tx_discontinuation1", label = "TP on to off treatment up to period 2", min = 0, max = 1, value = 0.069), 
      sliderInput(inputId = "tx_duration", label = "maximum treatment duration", min = 0, max = 50, value = 29), 
      numericInput(inputId = "c.Tx", label = "Cost of treatment", value = 26500 + (52/2)*78.35), 
      numericInput(inputId = "c.Tx_start", label = "Additional costs first cycle", value = 261.10*4 + 261.10*3*0.215), 
      numericInput(inputId = "wtp", label = "Willingness to pay threshold", min = 0, value = 100000), 
      numericInput(inputId = "discount_EFFECT", label = "Discount rate effects", min = 0, max = 1, value = 0.03), 
      numericInput(inputId = "discount_QALY", label = "Discount rate QALYs", min = 0, max = 1, value = 0.03), 
      numericInput(inputId = "discount_COST", label = "Discount rate costs", min = 0, max = 1, value = 0.03), 
      
      h3("Advanced inputs"), 
      
      numericInput(inputId = "p.mci_mod", label = "TP MCI to moderate dementia", min = 0, max = 1, value = 0),
      numericInput(inputId = "p.mci_sev", label = "TP MCI to severe dementia", min = 0, max = 1, value = 0),
      numericInput(inputId = "p.mil_mci", label = "TP mild dementia to MCI", min = 0, max = 1, value = 0.03),
      numericInput(inputId = "p.mil_sev", label = "TP mild dementia to severe dementia", min = 0, max = 1, value = 0.04),
      numericInput(inputId = "p.mod_mil", label = "TP moderate dementia to mild dementia", min = 0, max = 1, value = 0.03),
      numericInput(inputId = "p.mod_sev", label = "TP moderate dementia to severe dementia", min = 0, max = 1, value = 0.42),
      numericInput(inputId = "p.sev_mil", label = "TP severe dementia to mild dementia", min = 0, max = 1, value = 0),
      numericInput(inputId = "p.sev_mod", label = "TP severe dementia to moderate dementia", min = 0, max = 1, value = 0.02),
      numericInput(inputId = "p.mci_i", label = "TP MCI to institution care setting", min = 0, max = 1, value = 0.024),
      numericInput(inputId = "p.mil_i", label = "TP mild to institution care setting", min = 0, max = 1, value = 0.038),
      numericInput(inputId = "p.mod_i", label = "TP moderate to institution care setting", min = 0, max = 1, value = 0.110),
      numericInput(inputId = "p.sev_i", label = "TP severe to institution care setting", min = 0, max = 1, value = 0.259),
      numericInput(inputId = "hr.mort_mci", label = "HR death MCI compared to general population", min = 0, max = 1, value = 1.82),
      numericInput(inputId = "hr.mort_mil", label = "HR death mild dementia compared to general population", min = 0, max = 1, value = 2.92),
      numericInput(inputId = "hr.mort_mod", label = "HR death moderate dementia compared to general population", min = 0, max = 1, value = 3.85),
      numericInput(inputId = "hr.mort_sev", label = "HR death severe dementia compared to general population", min = 0, max = 1, value = 9.52),
      numericInput(inputId = "rr.tx_mci_mod", label = "RR treatment effect MCI to moderate dementia", min = 0, max = 1, value = 1),
      numericInput(inputId = "rr.tx_mci_sev", label = "RR treatment effect MCI to severe dementia", min = 0, max = 1, value = 1),
      numericInput(inputId = "rr.tx_mci_mil_dis", label = "RR treatment effect for transition MCI to mild dementia when no longer on treatment", min = 0, max = 1, value = 1), 
      numericInput(inputId = "rr.tx_mci_mod_dis", label = "RR treatment effect for transition MCI to moderate dementia when no longer on treatment", min = 0, max = 1, value = 1), 
      numericInput(inputId = "rr.tx_mci_sev_dis", label = "RR treatment effect for transition MCI to severe dementia when no longer on treatment", min = 0, max = 1, value = 1), 
      numericInput(inputId = "rr.tx_mil_mod_dis", label = "RR treatment effect for transition mild dementia to moderate dementia when no longer on treatment", min = 0, max = 1, value = 1), 
      numericInput(inputId = "rr.tx_mil_sev_dis", label = "RR treatment effect for transition mild dementia to severe dementia when no longer on treatment", min = 0, max = 1, value = 1), 
      numericInput(inputId = "p.tx_discontinuation2", label = "TP from on to off treatment from period 2 up to maximum treatment duration", min = 0, max = 1, value = 0), 
      numericInput(inputId = "tx_discontinuation2_begin", label = "Cycle number at which discontinuation period 2 starts", min = 0, max = 50, value = 2), 
      numericInput(inputId = "tx_waning", label = "Treatment effect waning per cycle expressed as relative reduction in treatment effect when on treatment", min = 0, max = 1, value = 0), 
      numericInput(inputId = "tx_waning_dis", label = "Treatment effect waning per cycle expressed as relative reduction in treatment effect when no longer on treatment", min = 0, max = 1, value = 0), 
      numericInput(inputId = "u.mci_pt", label = "Utility patient in state MCI community setting", value = 0.851 - 0.17),
      numericInput(inputId = "u.mil_pt", label = "Utility patient in state mild dementia community setting", value = 0.851 - 0.22),
      numericInput(inputId = "u.mod_pt", label = "Utility patient in state moderate dementia community setting", value = 0.851 - 0.36),
      numericInput(inputId = "u.sev_pt", label = "Utility patient in state severe dementia community setting", value = 0.851 - 0.53),
      numericInput(inputId = "u.mci_pt_i", label = "Utility patient in state MCI institution setting", value = 0.851 - 0.17),
      numericInput(inputId = "u.mil_pt_i", label = "Utility patient in state mild dementia institution setting", value = 0.851 - 0.19),
      numericInput(inputId = "u.mod_pt_i", label = "Utility patient in state moderate dementia institution setting", value = 0.851 - 0.42),
      numericInput(inputId = "u.sev_pt_i", label = "Utility patient in state severe dementia institution setting", value = 0.851 - 0.59),
      numericInput(inputId = "u.mci_ic", label = "Utility informal caregiver in state MCI community setting", value = -0.03),
      numericInput(inputId = "u.mil_ic", label = "Utility informal caregiver in state mild dementia community setting", value = -0.05),
      numericInput(inputId = "u.mod_ic", label = "Utility informal caregiver in state moderate dementia community setting", value = -0.08),
      numericInput(inputId = "u.sev_ic", label = "Utility informal caregiver in state severe dementia community setting", value = -0.10),
      numericInput(inputId = "u.mci_ic_i", label = "Utility informal caregiver in state MCI institution setting", value = -0.03),
      numericInput(inputId = "u.mil_ic_i", label = "Utility informal caregiver in state mild dementia institution setting", value = -0.05),
      numericInput(inputId = "u.mod_ic_i", label = "Utility informal caregiver in state moderate dementia institution setting", value = -0.08),
      numericInput(inputId = "u.sev_ic_i", label = "Utility informal caregiver in state severe dementia institution setting", value = -0.10),
      numericInput(inputId = "u.Tx_start", label = "Additional utility patient in first cycle (not half-cycle corrected)", value = -0.14 * (12/52) * 0.035),
      numericInput(inputId = "c.mci_hc", label = "Cost patient health care in state mci community setting", value = 6042*1.12 +  460),
      numericInput(inputId = "c.mil_hc", label = "Cost patient health care in state mild dementia community setting", value = 6042*1.56 +  965 + 0.21*365*0.333),
      numericInput(inputId = "c.mod_hc", label = "Cost patient health care in state moderate dementia community setting", value = 6042*1.93 + 1544 + 0.66*365*0.333),
      numericInput(inputId = "c.sev_hc", label = "Cost patient health care in state severe dementia community setting", value = 6042*1.93 + 1930),
      numericInput(inputId = "c.mci_hc_i", label = "Cost informal care in state mci institution setting", value = 6042*1.12 +  460),
      numericInput(inputId = "c.mil_hc_i", label = "Cost informal care in state mild dementia institution setting", value = 6042*1.56 +  965 + 0.21*365*0.333),
      numericInput(inputId = "c.mod_hc_i", label = "Cost informal care in state moderate dementia institution setting", value = 6042*1.93 + 1544 + 0.66*365*0.333),
      numericInput(inputId = "c.sev_hc_i", label = "Cost informal care in state severe dementia institution setting", value = 6042*1.93 + 1930),
      numericInput(inputId = "c.mci_sc", label = "Costs patient social care in state mci community setting", value = 0),
      numericInput(inputId = "c.mil_sc", label = "Costs patient social care in state mild dementia community setting", value = 0),
      numericInput(inputId = "c.mod_sc", label = "Costs patient social care in state moderate dementia community setting", value = 0),
      numericInput(inputId = "c.sev_sc", label = "Costs patient social care in state severe dementia community setting", value = 0),
      numericInput(inputId = "c.mci_sc_i", label = "Costs patient social care in state mci institution setting", value = 7394*12),
      numericInput(inputId = "c.mil_sc_i", label = "Costs patient social care in state mild dementia institution setting", value = 7394*12),
      numericInput(inputId = "c.mod_sc_i", label = "Costs patient social care in state moderate dementia institution setting", value = 7394*12),
      numericInput(inputId = "c.sev_sc_i", label = "Costs patient social care in state severe dementia institution setting", value = 7394*12),
      numericInput(inputId = "c.mci_ic", label = "Cost informal care in state mci community setting", value =  69*12*32.46 + 0.204*0.049*20*52*32.46),
      numericInput(inputId = "c.mil_ic", label = "Cost informal care in state mild dementia community setting", value = 113*12*32.46 + 0.112*0.086*20*52*32.46),
      numericInput(inputId = "c.mod_ic", label = "Cost informal care in state moderate dementia community setting", value = 169*12*32.46),
      numericInput(inputId = "c.sev_ic", label = "Cost informal care in state severe dementia community setting", value = 298*12*32.46),
      numericInput(inputId = "c.mci_ic_i", label = "Cost informal care in state mci institution setting", value =  69*12*32.46*0.44 + 0.204*0.049*20*52*32.46),
      numericInput(inputId = "c.mil_ic_i", label = "Cost informal care in state mild dementia institution setting", value = 113*12*32.46*0.44 + 0.112*0.086*20*52*32.46),
      numericInput(inputId = "c.mod_ic_i", label = "Cost informal care in state moderate dementia institution setting", value = 169*12*32.46*0.44),
      numericInput(inputId = "c.sev_ic_i", label = "Cost informal care in state severe dementia institution setting", value = 298*12*32.46*0.44),
      checkboxInput(inputId = "half_cycle_correction", label = "Apply half-cycle correction", value = TRUE)
    ),
    
    
    
    ######################################## * MAIN PANEL ########################################
    
    mainPanel(
      p("This is the beta version of the IPECAD Open-Source Model (version 2) for cost-effectiveness analysis of Alzheimer’s disease interventions. For important background information see: "), 
      a("github.com/ronhandels/ipecad", href="http://github.com/ronhandels/ipecad"),
      helpText(
        ("Abbreviations:"), br(), 
        ("dth = death"), br(),
        ("int = intervention strategy"), br(),         
        ("LY = life years"), br(),         
        ("mci = mild cognitive impairment"), br(),         
        ("mil = mild dementia"), br(),         
        ("mod = moderate dementia"), br(),         
        ("NHB = net health benefit"), br(),         
        ("QALY = quality-adjusted life years"), br(),         
        ("sev = severe dementia"), br(), 
        ("soc = standard of care strategy")
      ), 
      
      #tableOutput("test"), # for testing purposes
      
      h2("Summary"),
      tableOutput("table.summary"),
      
      h2("Mean time in state"),
      plotOutput("plot.timestate"),
      
      h2("State trace (table):"), 
      h3("Standard of care strategy"), 
      tableOutput("table.tracesoc"), 
      h3("Intervention strategy"), 
      tableOutput("table.traceint"), 
      
      h2("State trace (plot): plot"),
      plotOutput("plot.trace"),
      
      h2("Incremental cost-effectiveness ratio"),
      tableOutput("table.icer"),
      helpText("Effect=QALY; Inc=incremental"),

      h2("Incremental cost-effectiveness plane"),
      plotOutput("plot.icer")
      
    )
  )
)



######################################## SERVER ########################################

# Define server logic
server <- function(input, output, session) {
  
  # load libraries
  library(dampack) # load package
  
  
  
  ######################################## * RESULT ########################################
  
  # Create a reactive expression (inputs are called through 'input$...')
  result <- reactive({
    
    
    
    ######################################## ** PREPARE MODEL ########################################
    
    # put UI inputs into list of model inputs
    m.lifetable_US_2019 <- as.matrix(read.csv(file="lifetable_US_2019_ssa.csv", header=TRUE))[,c("male","female")] # load life table
    m.mortality_rate_US_2019 <- -log(1-(m.lifetable_US_2019)) # convert probability to rate
    ## weight rate for male and female
    m.mortality_rate_US_2019 <- cbind(m.mortality_rate_US_2019, weighted=NA)
    m.mortality_rate_US_2019[,"weighted"] <- m.mortality_rate_US_2019[,"male"] * 0.48 + m.mortality_rate_US_2019[,"female"] * 0.52
    
    l.inputs <- list(
      v.names_state = c("mcion_c","mciof_c","milon_c","milof_c","mod_c","sev_c","mci_i","mil_i","mod_i","sev_i","dth"), # disease states: mci = mild cognitive impairment; mil = mild dementia; mod = moderate dementia; sev = severe dementia; dth = dead; x_i = living in institutional setting (without '_i' = living in community)
      v.names_strat = c("soc","int"), 
      age_start = input$age_start, 
      sex = input$sex, 
      p.starting_state_mci = input$p.starting_state_mci, 
      n.cycle = input$n.cycle, 
      p.mci_mil = input$p.mci_mil, 
      p.mci_mod = input$p.mci_mod, 
      p.mci_sev = input$p.mci_sev, 
      p.mil_mci = input$p.mil_mci, 
      p.mil_mod = input$p.mil_mod, 
      p.mil_sev = input$p.mil_sev, 
      p.mod_mil = input$p.mod_mil, 
      p.mod_sev = input$p.mod_sev, 
      p.sev_mil = input$p.sev_mil, 
      p.sev_mod = input$p.sev_mod, 
      p.mci_i = input$p.mci_i, 
      p.mil_i = input$p.mil_i, 
      p.mod_i = input$p.mod_i, 
      p.sev_i = input$p.sev_i, 
      m.r.mortality = m.mortality_rate_US_2019, 
      hr.mort_mci = input$hr.mort_mci, 
      hr.mort_mil = input$hr.mort_mil, 
      hr.mort_mod = input$hr.mort_mod, 
      hr.mort_sev = input$hr.mort_sev, 
      rr.tx_mci_mil = input$rr.tx_mci_mil, 
      rr.tx_mci_mod = input$rr.tx_mci_mod, 
      rr.tx_mci_sev = input$rr.tx_mci_sev, 
      rr.tx_mil_mod = input$rr.tx_mil_mod, 
      rr.tx_mil_sev = input$rr.tx_mil_sev, 
      rr.tx_mci_mil_dis = input$rr.tx_mci_mil_dis, 
      rr.tx_mci_mod_dis = input$rr.tx_mci_mod_dis, 
      rr.tx_mci_sev_dis = input$rr.tx_mci_sev_dis, 
      rr.tx_mil_mod_dis = input$rr.tx_mil_mod_dis, 
      rr.tx_mil_sev_dis = input$rr.tx_mil_sev_dis, 
      p.tx_discontinuation1 = input$p.tx_discontinuation1, 
      p.tx_discontinuation2 = input$p.tx_discontinuation2, 
      tx_discontinuation2_begin = input$tx_discontinuation2_begin, 
      tx_duration = input$tx_duration, 
      tx_waning = input$tx_waning, 
      tx_waning_dis = input$tx_waning_dis, 
      u.mci_pt = input$u.mci_pt, 
      u.mil_pt = input$u.mil_pt, 
      u.mod_pt = input$u.mod_pt, 
      u.sev_pt = input$u.sev_pt, 
      u.mci_pt_i = input$u.mci_pt_i, 
      u.mil_pt_i = input$u.mil_pt_i, 
      u.mod_pt_i = input$u.mod_pt_i, 
      u.sev_pt_i = input$u.sev_pt_i, 
      u.mci_ic = input$u.mci_ic, 
      u.mil_ic = input$u.mil_ic, 
      u.mod_ic = input$u.mod_ic, 
      u.sev_ic = input$u.sev_ic, 
      u.mci_ic_i = input$u.mci_ic_i, 
      u.mil_ic_i = input$u.mil_ic_i, 
      u.mod_ic_i = input$u.mod_ic_i, 
      u.sev_ic_i = input$u.sev_ic_i, 
      u.Tx_start = input$u.Tx_start, 
      c.mci_hc = input$c.mci_hc, 
      c.mil_hc = input$c.mil_hc, 
      c.mod_hc = input$c.mod_hc, 
      c.sev_hc = input$c.sev_hc, 
      c.mci_hc_i = input$c.mci_hc_i, 
      c.mil_hc_i = input$c.mil_hc_i, 
      c.mod_hc_i = input$c.mod_hc_i, 
      c.sev_hc_i = input$c.sev_hc_i, 
      c.mci_sc = input$c.mci_sc, 
      c.mil_sc = input$c.mil_sc, 
      c.mod_sc = input$c.mod_sc, 
      c.sev_sc = input$c.sev_sc, 
      c.mci_sc_i = input$c.mci_sc_i, 
      c.mil_sc_i = input$c.mil_sc_i, 
      c.mod_sc_i = input$c.mod_sc_i, 
      c.sev_sc_i = input$c.sev_sc_i, 
      c.mci_ic =  input$c.mci_ic, 
      c.mil_ic = input$c.mil_ic, 
      c.mod_ic = input$c.mod_ic, 
      c.sev_ic = input$c.sev_ic, 
      c.mci_ic_i =  input$c.mci_ic_i, 
      c.mil_ic_i = input$c.mil_ic_i, 
      c.mod_ic_i = input$c.mod_ic_i, 
      c.sev_ic_i = input$c.sev_ic_i, 
      c.Tx = input$c.Tx, 
      c.Tx_start = input$c.Tx_start, 
      discount_EFFECT = input$discount_EFFECT, 
      discount_QALY = input$discount_QALY, 
      discount_COST = input$discount_COST, 
      wtp = input$wtp, 
      half_cycle_correction = input$half_cycle_correction
    )
    
    # load model functions
    source("functions_model.R")
    
    
    
    ######################################## ** RUN MODEL ########################################
    
    # run scenario and results
    l.out <- f.run_scenario(l.inputs = l.inputs, detailed = TRUE)
    
    # additional results
    m.result <- matrix(data = NA, nrow = ncol(l.out$l.out$soc$m.out), ncol = 3, dimnames = list(colnames(l.out$l.out$soc$m.out), c("soc","int","dif")))
    m.result[,"soc"] <- colSums(l.out$l.out$soc$m.out)
    m.result[,"int"] <- colSums(l.out$l.out$int$m.out)
    m.result[,"dif"] <- m.result[,"int"] - m.result[,"soc"]
    
    
    
    ######################################## ** PREPARE OUTCOMES ########################################
    
    # test (for testing purposes)
    test <- input$p.mci_mil
    
    # table: summary
    table.summary_data <- m.result[c("mci","mil","mod","sev","ly","qaly","cost"),c("soc","int","dif")]
    table.summary <- format(
      table.summary_data, 
      digits = 2, 
      scientific = FALSE, 
      big.mark = ","
    )
    
    # plot: time in state
    plot.timestate_data <- m.result[c("mci","mil","mod","sev"),c("int","soc")]
    
    # table, plot: state trace SOC
    tableplot.tracesoc_data <- cbind(
      mci = rowSums(l.out[["l.out"]][["soc"]][["m.trace"]][,c("mcion_c","mciof_c","mci_i")]), 
      mil = rowSums(l.out[["l.out"]][["soc"]][["m.trace"]][,c("milon_c","milof_c","mil_i")]),
      mod = rowSums(l.out[["l.out"]][["soc"]][["m.trace"]][,c("mod_c","mod_i")]), 
      sev = rowSums(l.out[["l.out"]][["soc"]][["m.trace"]][,c("sev_c","sev_i")]), 
      dth =         l.out[["l.out"]][["soc"]][["m.trace"]][,c("dth")]
    )
    
    # table, plot: state trace INT
    tableplot.traceint_data <- cbind(
      mci = rowSums(l.out[["l.out"]][["int"]][["m.trace"]][,c("mcion_c","mciof_c","mci_i")]), 
      mil = rowSums(l.out[["l.out"]][["int"]][["m.trace"]][,c("milon_c","milof_c","mil_i")]),
      mod = rowSums(l.out[["l.out"]][["int"]][["m.trace"]][,c("mod_c","mod_i")]), 
      sev = rowSums(l.out[["l.out"]][["int"]][["m.trace"]][,c("sev_c","sev_i")]), 
      dth =         l.out[["l.out"]][["int"]][["m.trace"]][,c("dth")]
    )
    
    # table, plot: icer
    icer <- calculate_icers(
      cost = l.out[["df.out"]][,"COST"],
      effect = l.out[["df.out"]][,"QALY"],
      strategies = l.out[["df.out"]][,"strategy"]
    )
    plot.icer <- icer
    table.icer <- format(
      icer,
      digits = 2,
      scientific = FALSE,
      big.mark = ","
    )
    
    # plot: cost difference by sector over time
    m.cost_incr_pos <- m.cost_incr_neg <- l.out[["l.out"]][["int"]][["m.out"]][,c("cost_dx","cost_tx","cost_hc","cost_sc","cost_ic")] - l.out[["l.out"]][["soc"]][["m.out"]][,c("cost_dx","cost_tx","cost_hc","cost_sc","cost_ic")] # split positive and negative
    m.cost_incr_pos[m.cost_incr_pos<0] <- 0
    m.cost_incr_neg[m.cost_incr_neg>=0] <- 0
    plot.incr_pos <- m.cost_incr_pos
    plot.incr_neg <- m.cost_incr_neg
    
    
    
    ######################################## ** LIST OUTCOMES ########################################
    
    # create list of output for each table and plot
    list(
      test = test, 
      table.summary = table.summary, 
      plot.timestate_data = plot.timestate_data, 
      tableplot.tracesoc_data = tableplot.tracesoc_data, 
      tableplot.traceint_data = tableplot.traceint_data, 
      plot.icer = plot.icer,
      table.icer = table.icer, 
      plot.incr_pos = plot.incr_pos, 
      plot.incr_neg = plot.incr_neg
    )
    
  })
  
  
  
  ######################################## * RENDER OUTCOMES ########################################
  
  # render test
  output$test <- renderTable(
    {
      result()[["test"]]
    }
  )
  
  # table: summary
  output$table.summary <- renderTable(
    {
      result()[["table.summary"]]
    },
    rownames = TRUE
  )
  
  # plot: time in state
  output$plot.timestate <- renderPlot(
    {
      plot.timestate_data <- result()[["plot.timestate_data"]]
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
    }
  )
  
  # table: state trace SOC
  output$table.tracesoc <- renderTable(
    {
      tableplot.tracesoc_data <- result()[["tableplot.tracesoc_data"]]
      round(tableplot.tracesoc_data[1:min(nrow(tableplot.tracesoc_data),10),],2)
    }
  )
  
  # table: state trace INT
  output$table.traceint <- renderTable(
    {
      tableplot.traceint_data <- result()[["tableplot.traceint_data"]]
      round(tableplot.traceint_data[1:min(nrow(tableplot.traceint_data),10),],2)
    }
  )
  
  # plot: state trace
  output$plot.trace <- renderPlot(
    {
      tableplot.tracesoc_data <- result()[["tableplot.tracesoc_data"]]
      tableplot.traceint_data <- result()[["tableplot.traceint_data"]]
      
      xrange <- 0:(nrow(tableplot.tracesoc_data)-1)
      xx <- c(xrange, rev(xrange)) # prepare polygon x-values
      yy_mci <- c(tableplot.tracesoc_data[,"mci"], rev(tableplot.traceint_data[,"mci"])) # polygon y-values
      yy_mil <- c(tableplot.tracesoc_data[,"mil"], rev(tableplot.traceint_data[,"mil"])) # idem
      yy_mod <- c(tableplot.tracesoc_data[,"mod"], rev(tableplot.traceint_data[,"mod"])) # idem
      yy_sev <- c(tableplot.tracesoc_data[,"sev"], rev(tableplot.traceint_data[,"sev"])) # idem
      yy_dth <- c(tableplot.tracesoc_data[,"dth"], rev(tableplot.traceint_data[,"dth"])) # idem
      
      par(mar=c(5, 4, 4, 1)+0.1, xpd=FALSE)
      matplot(
        x = xrange, 
        y = cbind(tableplot.tracesoc_data,tableplot.traceint_data), 
        type = "n", 
        xlab = "cycle", 
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
        x = xrange, 
        y = tableplot.tracesoc_data, 
        type = "l",
        lty = 1,
        col = c("green","yellow","orange","red","black")
      )
      matlines(
        x = xrange, 
        y = tableplot.traceint_data, 
        type = "l",
        lty = 2,
        col = c("green","yellow","orange","red","black")
      )
      legend(x="topright", legend=c("mci","mil","mod","sev","dth"), col=c("green","yellow","orange","red","black"), lty=1, bg="white")
      legend(x="right", legend=c("soc","int"), col="black", lty=c(1,2), bg="white")
    }
  )
  
  # table: icer
  output$table.icer <- renderTable(
    {
      result()[["table.icer"]]
    },
    rownames = TRUE
  )

  # plot: icer
  output$plot.icer <- renderPlot(
    {
      par(mar=c(5, 4, 4, 2)+0.1, xpd=FALSE)
      plot(result()[["plot.icer"]], label="all")
    }
  )
  
  # plot: cost difference by sector over time
  output$plot.icer <- renderPlot(
    {
      plot.incr_pos <- result()[["plot.incr_pos"]]
      plot.incr_neg <- result()[["plot.incr_neg"]]
      
      par(mar=c(5, 4, 4, 2)+0.1, xpd=FALSE)
      barplot(
        height = t(plot.incr_pos),
        beside = F,
        xlab = "cycle",
        ylab = "incremental costs",
        ylim = c(
          min(plot.incr_neg) + min(plot.incr_neg)*0.10, 
          max(plot.incr_pos) + max(plot.incr_pos)*0.10
        ),
        col = rainbow(5), 
        names.arg = 1:nrow(plot.incr_pos),
        main = "costs by sector over time"
      )
      barplot(
        height = t(plot.incr_neg),
        beside = F,
        col = rainbow(5),
        add = T
      )
      legend(x = "topright", legend = c("diagnostic","treatment","health","social","informal"), fill = rainbow(5))
    }
  )
  
}

# Create Shiny app
shinyApp(ui, server)