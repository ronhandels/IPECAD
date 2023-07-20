library(shiny)

######################################## UI ########################################
ui <- fluidPage(
  
  # App title
  titlePanel("open-source model"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    
    
    ######################################## * SIDEBAR PANEL ########################################
    sidebarPanel(
      p("THIS IS AN ALPHA VERSION FOR TESTING ONLY. For details, see "), 
      a("github.com/ronhandels/ipecad", href="http://github.com/ronhandels/ipecad"),
      p(" "), 
      
      sliderInput(inputId="in_age_start_end", label="age start and end:", min = 50, max = 99, value = c(70,99)), 
      radioButtons(inputId="sex", label="sex:", choices=c("male","female")), 
      sliderInput(inputId="p.mci_mil", label="mci>mil (transition probability):", min = 0.01, max = 0.90, value = 0.206), 
      sliderInput(inputId="rr.tx", label="treatment effect (RR):", min = 0.00, max = 1.00, value = 0.75), 
      sliderInput(inputId="p.discontinuation", label="discontinuation (proportion per year):", min = 0.00, max = 1.00, value = 0.10), 
      sliderInput(inputId="tx_duration", label="treatment duration (years):", min = 0, max = 20, value = 7), 
      sliderInput(inputId="p.starting_state_mci", label="starting state MCI (proportion):", min = 0.00, max = 1.00, value = 1.00), 
      sliderInput(inputId="c.Tx", label="costs treatment (USD):", min=0, max=50000, value=5000, step=1000), 
      # selectInput(inputId="p.mci_mil", label="TP mci>mil:", choices=c(default_0.206=0.206, vos1_0.15=0.15), selectize=TRUE), 
      # sliderInput(inputId="in_age_start", label="starting age:", min=50, max=99, value=70), 
      # sliderInput(inputId="decimal", label="Decimal:", min = 0, max = 1, value = 0.5, step = 0.1), 
      # sliderInput(inputId="format", label="Custom Format:", min = 0, max = 10000, value = 0, step = 2500, pre = "$", sep = ",", animate = TRUE), # custom currency format for with basic animation
      # sliderInput(inputId="animation", label="Looping Animation:", min = 1, max = 2000, value = 1, step = 10, animate = animationOptions(interval = 300, loop = TRUE)), # Input: Animation with custom interval (in ms) to control speed, plus looping
      helpText("abbreviations: dif=difference; dth=death; int=intervention; mci=mild cognitive impairment; mil=mild dementia; RR=relative risk; soc=standard of care; TP=transition probability; ")
    ),
    
    
    ######################################## * MAIN PANEL ########################################
    mainPanel(
      h1("Model outcomes"),
      p("THIS IS AN ALPHA VERSION FOR TESTING ONLY. For details, see "), 
      a("github.com/ronhandels/ipecad", href="http://github.com/ronhandels/ipecad"),
      p(" "), 
      
      tableOutput("table_summary"),
      plotOutput("trace"),
      plotOutput("plot2", width = "100%"),
      tableOutput("table_icer"), 
      plotOutput("icer")
      
    )
  )
)



######################################## SERVER ########################################

# Define server logic for slider examples
server <- function(input, output, session) {

  # load libraries
  library(dampack) # load package
  
  
  ######################################## * RESULT ########################################
  
  # Create a reactive expression (inputs are called through 'inputs$...')
  result <- reactive({
    
    
    ######################################## ** PREPARE MODEL ########################################
    
    # put UI inputs into list of model inputes
    m.lifetable_US <- as.matrix(read.csv(file="lifetable_US.csv", header=TRUE))[,c("male","female")] # import life table and select only men/women and drop age column (make sure age corresponds to row number, i.e., start with age = 1)
    m.mortality_rate_US <- -log(1-(m.lifetable_US)) # convert probability to rate
    l.inputs <- list(
      v.names_state = c("mcion","mciof","milon","milof","mod","sev","dth"), # disease states: mci = mild cognitive impairment; mil = mild dementia; mod = moderate dementia; sev = severe dementia; dth = dead
      v.names_strat = c("soc","int"), # strategies: soc = standard of care strategy; int = intervention strategy
      age_start = input$in_age_start_end[[1]], # start age (dicated by the benchmark scenario)
      age_end = input$in_age_start_end[[2]], # end age (set at max)
      sex = input$sex, # sex of starting population (dependent is mortality table)
      p.mci_mil = as.numeric(input$p.mci_mil), # p.x_x: transition probability between states [Wimo, 2020: https://doi.org/10.3233/jad-191055]
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
      rr.tx_mci_mil = input$rr.tx, # treatment effect expressed as hazard ratio on transition rate
      rr.tx_mci_mod = 1, # idem
      rr.tx_mci_sev = 1, # idem
      rr.tx_mil_mod = input$rr.tx, # assumed same effect as for MCI to dementia
      tx_waning = 0.05, # assumed annual waning of treatment
      p.discontinuation1 = input$p.discontinuation, # discontinuation at year 1
      p.discontinuation_x = input$p.discontinuation, # annual proportion discontinuation
      tx_duration = input$tx_duration, # maximum treatment duration
      p.starting_state_mci = input$p.starting_state_mci, # proportion starting in disease state MCI, remaining from 1 will start in 'mil' (all will start as 'of' in 'soc' and 'on' in 'int')
      u.mci = 0.73, # u.x: utility in state [https://doi.org/10.1016/j.jalz.2019.05.004]
      u.mil = 0.69, # idem
      u.mod = 0.53, # idem
      u.sev = 0.38, # idem
      c.mci = (1254 +  222) * 12 * (1-0    ) + (1254 + 8762) * 12 * 0, # c.x: costs in state, build up as montly costs in patient health and social care by care setting (community/residential) multiplied by 12 (annual costs) [Tahami, 2023: https://doi.org/10.1007/s40120-023-00460-1 table 2] and multiplied by proportion in setting [Tahami, 2023: https://doi.org/10.1007/s40120-023-00460-1 table 1]
      c.mil = (1471 +  410) * 12 * (1-0.038) + (1471 + 8762) * 12 * 0.038, # idem
      c.mod = (1958 +  653) * 12 * (1-0.110) + (1958 + 8762) * 12 * 0.110, # idem
      c.sev = (2250 + 1095) * 12 * (1-0.259) + (2250 + 8762) * 12 * 0.259, # idem
      c.Tx = input$c.Tx, # treatment costs
      c.Tx_diagnostics1 = 2000, # costs diagnostics cycle 1 (not half-cycle corrected)
      discount_QALY = 0.035, # discount rate
      discount_COST = 0.035, # # discount rate
      wtp = 40000 # willingness to pay
    )
    
    # load model functions
    source("functions_model.R")
    
    
    ######################################## ** RUN MODEL ########################################
    out_base <- f.run_scenario(l.inputs = l.inputs, detailed = TRUE)
    
    
    ######################################## ** PREPARE OUTCOMES ########################################
    
    # test (not used)
    test <- input$p.mci_mil
    
    ## table: summary outcomes
    temp_table1 <- rbind(
      out_base[["df.out_sum"]][,-1], 
      incremental=out_base[["df.out_sum"]]["int",-1] - out_base[["df.out_sum"]]["soc",-1] # calculate difference between 'soc' and 'int' strategies
    )
    df.table1 <- temp_table1
    df.summary <- format(temp_table1, digits=2, scientific=FALSE, big.mark=",") # format 
    
    ## plot: state trace
    m.trace_soc <- cbind(
      mci=out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"mcion"] + out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"mciof"], 
      mil=out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"milon"] + out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"milof"], 
      mod=out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"mod"], 
      sev=out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"sev"], 
      dth=out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"dth"]
    )
    m.trace_int <- cbind(
      mci=out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"mcion"] + out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"mciof"], 
      mil=out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"milon"] + out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"milof"], 
      mod=out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"mod"], 
      sev=out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"sev"], 
      dth=out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"dth"]
    )
    a.trace <- array(data=c(m.trace_soc,m.trace_int), dim=c(nrow(m.trace_soc),5,2), dimnames=list(NULL,colnames(m.trace_soc),c("soc","int")))
    
    ## plot: mean time in state
    m.plot2 <- cbind(soc=colSums(m.trace_soc), int=colSums(m.trace_int))
    m.plot2 <- m.plot2[c("mci","mil","mod","sev"),]
    
    ## plot: ICE-plane
    icer <- calculate_icers(
      cost = out_base[["df.out_sum"]][,"COST"],
      effect = out_base[["df.out_sum"]][,"QALY"],
      strategies = out_base[["df.out_sum"]][,"strategy"]
    )
    
    ## table: icer
    df.icer <- icer
    df.icer <- format(df.icer, digits=2, scientific=FALSE, big.mark=",")

    
    ######################################## ** LIST OUTCOMES ########################################
    
    # create list of output for each of the tables/figures
    list(
      test = test,
      df.summary = df.summary, 
      a.trace = a.trace, 
      plot2 = m.plot2,
      icer = icer, 
      df.icer = df.icer
    )
    
  })
  
  
  
  ######################################## * RENDER OUTCOMES ########################################
  
  # # render test
  # output$test <- renderTable({
  #   result()[["test"]]
  # })
  
  # table: summary outcomes
  output$table_summary <- renderTable(
    { 
      result()[["df.summary"]]
    }, 
    rownames = TRUE
  )
  
  # plot: state trace
  output$trace <- renderPlot(
    {
      a.trace <- result()[["a.trace"]]
      v.age_range <- c(input$in_age_start_end[[1]]:(input$in_age_start_end[[2]]-1)) # store age range
      xx <- c(v.age_range, rev(v.age_range)) # prepare polygon x-values
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
        y = a.trace[,,"soc"], 
        type = "l",
        lty = 1,
        col = c("green","yellow","orange","red","black")
      )
      matlines(
        x = v.age_range, 
        y = a.trace[,,"int"], 
        type = "l",
        lty = 2,
        col = c("green","yellow","orange","red","black")
      )
      legend(x="topright", legend=c("mci","mil","mod","sev","dth"), col=c("green","yellow","orange","red","black"), lty=1)
      legend(x="right", legend=c("soc","int"), col="black", lty=c(1,2))
    }
  )
  
  # plot: mean time in state
  output$plot2 <- renderPlot(
    {
      m.plot2 <- result()[["plot2"]]
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
    }
  )
  
  # plot: icer
  output$plot3 <- renderPlot(
    {
      par(mar=c(5, 4, 4, 2)+0.1, xpd=FALSE)
      plot(result()[["icer"]], label="all")
    }
  )
  
  # table: icer
  output$table_icer <- renderTable(
    {
      result()[["df.icer"]]
    }, 
    rownames = TRUE
  )
  
  
  
  
}

# Create Shiny app
shinyApp(ui, server)