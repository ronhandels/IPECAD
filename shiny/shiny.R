library(shiny)

######################################## UI ########################################
ui <- fluidPage(
  
  # App title
  titlePanel("IPECAD open-source model"),
  
  # Sidebar layout with input and output definitions
  sidebarLayout(
    
    
    ######################################## * SIDEBAR PANEL ########################################
    sidebarPanel(
      p(""), 
      sliderInput(inputId="in_age_start_end", label="Age start and end:", min = 50, max = 99, value = c(70,99)), 
      radioButtons(inputId="sex", label="Sex:", choices=c("female","male")), 
      sliderInput(inputId="p.mci_mil", label="Transition probability MCI to mild dementia:", min = 0.0, max = 1.0, value = 0.21), 
      helpText("amyloid positive & injury positive = 0.28; amyloid positive & injury undetermined = 0.25 [Vos, 2015]"), 
      sliderInput(inputId="p.discontinuation", label="Discontinuation (proportion per year):", min = 0.00, max = 1.00, value = 0.10), 
      sliderInput(inputId="rr.tx", label="Treatment effect (RR):", min = 0.00, max = 1.00, value = 0.75), 
      helpText("applied to transition MCI to mild dementia and to mild to moderate dementia"), 
      sliderInput(inputId="tx_waning", label="Treatment waning:", min = 0.00, max = 1.00, value = 0.05), 
      helpText("proportion per year treatment no longer effective"), 
      sliderInput(inputId="tx_duration", label="Treatment duration (years):", min = 0, max = 20, value = 7), 
      sliderInput(inputId="p.starting_state_mci", label="Starting state MCI (proportion):", min = 0.00, max = 1.00, value = 1.00), 
      helpText("1 minus this proportion starts in mild dementia"), 
      sliderInput(inputId="c.Tx", label="Costs treatment (USD):", min=0, max=50000, value=10000, step=1000), 
      helpText("when on treatment")
      # sliderInput(inputId="in_age_start", label="starting age:", min=50, max=99, value=70), 
      # sliderInput(inputId="decimal", label="Decimal:", min = 0, max = 1, value = 0.5, step = 0.1), 
      # sliderInput(inputId="format", label="Custom Format:", min = 0, max = 10000, value = 0, step = 2500, pre = "$", sep = ",", animate = TRUE), # custom currency format for with basic animation
      # sliderInput(inputId="animation", label="Looping Animation:", min = 1, max = 2000, value = 1, step = 10, animate = animationOptions(interval = 300, loop = TRUE)), # Input: Animation with custom interval (in ms) to control speed, plus looping
    ),
    
    
    ######################################## * MAIN PANEL ########################################
    mainPanel(
      p("This is the beta version of the IPECAD Open-Source Model v2 - Single-Domain for cost-effectiveness analysis of Alzheimer’s disease interventions. For important background information see: "), 
      a("github.com/ronhandels/ipecad", href="http://github.com/ronhandels/ipecad"),
      helpText("abbreviations: dth=death; int=intervention strategy; LY=life years; mci=mild cognitive impairment; mil=mild dementia; mod=moderate dementia; NHB=net health benefit; QALY=quality-adjusted life years; sev=severe dementia; soc=standard of care strategy"), 
      
      h2("Health-economic outcomes"),
      tableOutput("table_summary"),

      h2("State trace: plot"),
      plotOutput("trace"),

      h2("State trace:"), 
      h3("Standard of care strategy"), 
      tableOutput("trace_soc"), 
      h3("Intervention strategy"), 
      tableOutput("trace_int"), 

      h2("State trace: mean time in state"), 
      plotOutput("plot2", width = "100%"),
      
      h2("Incremental cost-effectiveness plane"), 
      plotOutput("icer"), 

      h2("Incremental cost-effectiveness ratio"), 
      tableOutput("table_icer"), 
      helpText("Effect=QALY; Inc=incremental")
      
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
    
    # put UI inputs into list of model inputs
    m.lifetable_US <- as.matrix(read.csv(file="life_tables/lifetable_US.csv", header=TRUE))[,c("male","female")] # import life table and select only men/women and drop age column (make sure age corresponds to row number, i.e., start with age = 1)
    m.mortality_rate_US <- -log(1-(m.lifetable_US)) # convert probability to rate
    l.inputs <- list(
      v.names_state = c("mcion","mciof","milon","milof","mod","sev","mci_i","mil_i","mod_i","sev_i","dth"), # disease states: mci = mild cognitive impairment; mil = mild dementia; mod = moderate dementia; sev = severe dementia; dth = dead; x_i = living in institutional setting (without '_i' = living in community)
      v.names_strat = c("soc","int"), # strategies: soc = standard of care strategy; int = intervention strategy
      age_start = input$in_age_start_end[[1]], 
      age_end = input$in_age_start_end[[2]], 
      sex = input$sex, 
      p.mci_mil = as.numeric(input$p.mci_mil), 
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
      p.discontinuation1 = input$p.discontinuation, 
      p.discontinuation_x = input$p.discontinuation, 
      m.r.mortality = m.mortality_rate_US, 
      hr.mort_mci = 1, 
      hr.mort_verymilddem = 1.82, 
      hr.mort_mil = 1.318, 
      hr.mort_mod = 2.419, 
      hr.mort_sev = 4.267, 
      rr.tx_mci_mil = input$rr.tx, 
      rr.tx_mci_mod = 1, 
      rr.tx_mci_sev = 1, 
      rr.tx_mil_mod = input$rr.tx, 
      tx_waning = input$tx_waning, 
      tx_duration = input$tx_duration, 
      p.starting_state_mci = input$p.starting_state_mci, 
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
      c.Tx = input$c.Tx, 
      c.Tx_diagnostics1 = 2000, 
      discount_QALY = 0.035, 
      discount_COST = 0.035, 
      wtp = 40000 
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
      mci=out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"mcion"] + out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"mciof"] + out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"mci_i"], 
      mil=out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"milon"] + out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"milof"] + out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"mil_i"], 
      mod=out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"mod"] + out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"mod_i"], 
      sev=out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"sev"] + out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"sev_i"], 
      dth=out_base[["l.out_strategy"]][["soc"]][["m.trace"]][,"dth"]
    )
    m.trace_int <- cbind(
      mci=out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"mcion"] + out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"mciof"] + out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"mci_i"], 
      mil=out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"milon"] + out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"milof"] + out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"mil_i"], 
      mod=out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"mod"] + out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"mod_i"], 
      sev=out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"sev"] + out_base[["l.out_strategy"]][["int"]][["m.trace"]][,"sev_i"], 
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
  
  
  # table: state trace SOC
  output$trace_soc <- renderTable(
    {
      round(result()[["a.trace"]][1:10,,"soc"], 2)
    } 
  )

  
  # table: state trace INT
  output$trace_int <- renderTable(
    {
      round(result()[["a.trace"]][1:10,,"int"], 2)
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
  output$icer <- renderPlot(
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