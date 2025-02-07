# Load Required Libraries
library(shiny)
library(ggplot2)
library(dplyr)
library(car)  
library(stats)
library(ggpubr)  
library(multcomp)
library(bslib)

# Helper Function to Determine Significance Labels
significance_labels <- function(p_value) {
  if (is.na(p_value)) {
    return("ns")
  } else if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else {
    return("ns")
  }
}

# Format the UI with bslib theme customization
ui <- fluidPage(
  # Define a custom Bootstrap theme using bslib
  theme = bs_theme(                     
    bootswatch = "slate"  # Bootswatch theme for a modern look
  ),
  
  titlePanel("Cell Growth Analysis Tool"),
  
  sidebarLayout(
    sidebarPanel(
      HTML('<img src="logo.png" width="100%" height="auto">'),
      downloadButton("downloadTemplate", "Download Template"),
      fileInput("data", "Upload Growth Data (CSV)", accept = ".csv"),
      actionButton("recommendTest", "Recommend Best Test"),
      textOutput("recommendedTest"),
      selectInput("statTest", "Select Statistical Test", 
                  choices = c(
                    "None" = "none",
                    "T-test" = "t.test", 
                    "ANOVA" = "anova", 
                    "Welch's ANOVA" = "welch.anova", 
                    "Kruskal-Wallis" = "kruskal",
                    "Pairwise Wilcoxon" = "wilcox"
                  )
      ),
      # Inputs for plot customization
      textInput("plotTitle", "Plot Title", value = "Cell Growth Over Time"),
      textInput("xLabel", "X-axis Label", value = "Time (Hours)"),
      textInput("yLabel", "Y-axis Label", value = "Mean Growth (Millions of Cells)"),
      actionButton("goButton", "Run Analysis"),
      downloadButton("downloadPlot", "Download Plot"),
      br(),
      br(),
      downloadButton("downloadResults", "Download Text Results"),
      br(),
      br(),
      tags$p("Created by Andy Ring"),
      tags$p("Version 1.0.1 | February, 7th 2025")
    ), 
    
    mainPanel(
      layout_columns(
        card(card_header("Growth Plot"),
             plotOutput("growthPlot")),
        card(card_header("Statistical Test Results"),
             tableOutput("statsTable")),
        card(card_header("Regression Equations"),
             tableOutput("regressionEquations")),  # Table for regression equations
        card(card_header("Growth Calculations"),
             tableOutput("growthMetrics")),  # Table for doubling time and growth rate
        col_widths = c(12, 4, 4, 4)
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  
  # Function to Create Downloadable Template
  createTemplate <- function() {
    data.frame(
      time = c(0, 1, 2, 3, 4, 5),
      growth = rep(NA, 6),  
      group = rep("Control", 6)  
    )
  }
  
  # Download Handler for Template
  output$downloadTemplate <- downloadHandler(
    filename = function() { "growth_data_template.csv" },
    content = function(file) { write.csv(createTemplate(), file, row.names = FALSE) }
  )
  
  # Reactive Expression to Read Uploaded Data
  data <- reactive({
    req(input$data)
    tryCatch(
      {
        df <- read.csv(input$data$datapath)
        # Validate Columns
        required_cols <- c("time", "growth", "group")
        if (!all(required_cols %in% colnames(df))) {
          stop("Missing required columns. Please ensure the CSV has 'time', 'growth', and 'group' columns.")
        }
        # Ensure correct data types
        df$time <- as.numeric(df$time)
        df$growth <- as.numeric(df$growth)
        df$group <- as.factor(df$group)
        return(df)
      },
      error = function(e) {
        showModal(modalDialog(
          title = "Data Upload Error",
          paste("Error:", e$message),
          easyClose = TRUE
        ))
        return(NULL)
      }
    )
  })
  
  # Recommend Best Statistical Test
  observeEvent(input$recommendTest, {
    growthData <- data()
    req(growthData)
    
    groups <- unique(growthData$group)
    recommendation <- ""
    
    if (length(groups) == 1) {
      bestTest <- "None"
      recommendation <- "No statistical test recommended because there is only one group present."
    } else if (length(groups) == 2) {
      normTest <- tryCatch(shapiro.test(growthData$growth)$p.value, error = function(e) NA)
      varTest <- tryCatch(leveneTest(growth ~ group, data = growthData)$`Pr(>F)`[1], error = function(e) NA)
      
      if (!is.na(normTest) && normTest > 0.05) {
        if (!is.na(varTest) && varTest > 0.05) {
          bestTest <- "T-test"
          recommendation <- "T-test because the data has 2 groups and is normally distributed with equal variance."
        } else {
          bestTest <- "Welch's T-test"
          recommendation <- "Welch's T-test because the data has 2 groups and is normally distributed but with unequal variance."
        }
      } else {
        bestTest <- "Wilcoxon Test"
        recommendation <- "Wilcoxon Test because the data is not normally distributed and has 2 groups."
      }
    } else {
      normTest <- tryCatch(shapiro.test(growthData$growth)$p.value, error = function(e) NA)
      varTest <- tryCatch(leveneTest(growth ~ group, data = growthData)$`Pr(>F)`[1], error = function(e) NA)
      
      if (!is.na(normTest) && normTest > 0.05) {
        if (!is.na(varTest) && varTest > 0.05) {
          bestTest <- "ANOVA"
          recommendation <- "ANOVA because there are more than 2 groups and data is normally distributed with equal variance."
        } else {
          bestTest <- "Welch's ANOVA"
          recommendation <- "Welch's ANOVA because there are more than 2 groups and data is normally distributed but with unequal variance."
        }
      } else {
        bestTest <- "Kruskal-Wallis Test"
        recommendation <- "Kruskal-Wallis Test because the data is not normally distributed and has more than 2 groups."
      }
    }
    
    output$recommendedTest <- renderText({ 
      paste("Recommended Test: ", bestTest, "\n", recommendation) 
    })
  })
  
  # Use eventReactive to run the analysis when "Run Analysis" is clicked.
  analysisResults <- eventReactive(input$goButton, {
    req(data())
    growthData <- data()
    growthData$group <- as.factor(growthData$group)
    
    # Summarize Data
    summarizedData <- growthData %>%
      group_by(time, group) %>%
      summarise(
        mean_growth = mean(growth, na.rm = TRUE),
        stderr = sd(growth, na.rm = TRUE) / sqrt(n()),
        .groups = 'drop'
      )
    
    # Generate Polynomial Regression Equations (for display)
    regressionEquations <- summarizedData %>%
      group_by(group) %>%
      do({
        if(n_distinct(.$time) >= 3) {
          model <- tryCatch(lm(mean_growth ~ poly(time, 2), data = .), error = function(e) NULL)
          if (!is.null(model)) {
            coefs <- coef(model)
            if(length(coefs) == 3) {
              data.frame(
                group = unique(.$group),
                equation = sprintf("y = %.3fx² + %.3fx + %.3f", coefs[3], coefs[2], coefs[1])
              )
            } else {
              data.frame(
                group = unique(.$group),
                equation = "Insufficient data for polynomial regression."
              )
            }
          } else {
            data.frame(
              group = unique(.$group),
              equation = "Regression failed."
            )
          }
        } else {
          data.frame(
            group = unique(.$group),
            equation = "Insufficient data for polynomial regression."
          )
        }
      })
    
    # Calculate Doubling Time and Growth Rate (Exponential Model)
    doublingGrowth <- summarizedData %>%
      group_by(group) %>%
      do({
        if(min(.$mean_growth, na.rm = TRUE) > 0 && n_distinct(.$time) >= 2) {
          expModel <- tryCatch(lm(log(mean_growth) ~ time, data = .), error = function(e) NULL)
          if (!is.null(expModel)) {
            slope <- coef(expModel)["time"]
            growth_rate <- slope
            doubling_time <- log(2) / slope
            data.frame(
              group = unique(.$group),
              growth_rate = growth_rate,
              doubling_time = doubling_time
            )
          } else {
            data.frame(
              group = unique(.$group),
              growth_rate = NA,
              doubling_time = NA
            )
          }
        } else {
          data.frame(
            group = unique(.$group),
            growth_rate = NA,
            doubling_time = NA
          )
        }
      })
    
    # Perform Selected Statistical Test (skip if "None" is selected)
    statsResult <- switch(input$statTest,
                          "none" = NULL,
                          "t.test" = {
                            if(length(levels(growthData$group)) == 2){
                              t.test(growth ~ group, data = growthData)
                            } else {
                              NULL
                            }
                          },
                          "anova" = {
                            if(length(levels(growthData$group)) > 1){
                              aov(growth ~ group, data = growthData)
                            } else {
                              NULL
                            }
                          },
                          "welch.anova" = {
                            if(length(levels(growthData$group)) > 1){
                              oneway.test(growth ~ group, data = growthData, var.equal = FALSE)
                            } else {
                              NULL
                            }
                          },
                          "kruskal" = {
                            if(length(levels(growthData$group)) > 1){
                              kruskal.test(growth ~ group, data = growthData)
                            } else {
                              NULL
                            }
                          },
                          "wilcox" = {
                            if(length(levels(growthData$group)) > 1){
                              pairwise.wilcox.test(growthData$growth, growthData$group, p.adjust.method = "bonferroni")
                            } else {
                              NULL
                            }
                          }
    )
    
    # Handle Pairwise Comparisons for Annotation
    pairwise_comparisons <- NULL
    if (input$statTest == "anova") {
      anova_model <- aov(growth ~ group, data = growthData)
      tukey_results <- TukeyHSD(anova_model)$group
      groups <- levels(growthData$group)
      pairwise_comparisons <- matrix(NA, nrow = length(groups), ncol = length(groups), 
                                     dimnames = list(groups, groups))
      for (comparison in rownames(tukey_results)) {
        groups_pair <- unlist(strsplit(comparison, "-"))
        if (length(groups_pair) == 2) {
          g1 <- groups_pair[1]
          g2 <- groups_pair[2]
          p_adj <- tukey_results[comparison, "p adj"]
          pairwise_comparisons[g1, g2] <- p_adj
          pairwise_comparisons[g2, g1] <- p_adj
        }
      }
    } else if (input$statTest == "kruskal") {
      pairwise_comparisons <- statsResult$p.value
      if(!is.matrix(pairwise_comparisons)){
        pairwise_comparisons <- as.matrix(pairwise_comparisons)
      }
    } else if (input$statTest == "wilcox") {
      pairwise_comparisons <- statsResult$p.value
      if(!is.matrix(pairwise_comparisons)){
        pairwise_comparisons <- as.matrix(pairwise_comparisons)
      }
    }
    
    # Generate Static Annotations Based on Pairwise Comparisons
    annotations <- data.frame()
    if (!is.null(pairwise_comparisons)) {
      groups <- levels(growthData$group)
      group_combos <- combn(groups, 2, simplify = FALSE)
      offset <- max(summarizedData$stderr, na.rm = TRUE) * 2  
      annotation_counter <- 0
      
      for (combo in group_combos) {
        group1 <- combo[1]
        group2 <- combo[2]
        if(all(c(group1, group2) %in% rownames(pairwise_comparisons))){
          p_val <- pairwise_comparisons[group1, group2]
        } else {
          p_val <- NA
        }
        if (!is.na(p_val)) {
          annotation_label <- significance_labels(p_val)
          annotation_counter <- annotation_counter + 1
          x_position <- mean(range(summarizedData$time))
          y_position <- max(summarizedData$mean_growth, na.rm = TRUE) + offset * annotation_counter
          annotations <- rbind(annotations, data.frame(
            x = x_position,
            y = y_position,
            label = paste0(group1, " vs ", group2, ": ", annotation_label)
          ))
        }
      }
    }
    
    # Return all computed results as a list
    list(
      summarizedData = summarizedData,
      regressionEquations = regressionEquations,
      doublingGrowth = doublingGrowth,
      statsResult = statsResult,
      annotations = annotations
    )
  })
  
  # Define Outputs using the results from analysisResults()
  output$growthPlot <- renderPlot({
    res <- analysisResults()
    req(res)
    
    p <- ggplot(res$summarizedData, aes(x = time, y = mean_growth, color = group)) +
      geom_point(size = 3) +  
      geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, linewidth = 1) +  
      geom_errorbar(aes(ymin = mean_growth - stderr, ymax = mean_growth + stderr), width = 0.2) +
      labs(
        title = input$plotTitle,
        x = input$xLabel,
        y = input$yLabel,
        color = "Group"
      ) +
      theme_light()
    
    if (nrow(res$annotations) > 0) {
      p <- p + geom_text(
        data = res$annotations, 
        aes(x = x, y = y, label = label), 
        size = 6, 
        vjust = -0.5,
        inherit.aes = FALSE
      )
    }
    
    p
  })
  
  output$regressionEquations <- renderTable({
    res <- analysisResults()
    req(res)
    res$regressionEquations
  })
  
  output$growthMetrics <- renderTable({
    res <- analysisResults()
    req(res)
    res$doublingGrowth
  })
  
  output$statsTable <- renderTable({
    # If "None" is selected, display a message.
    if (input$statTest == "none") {
      return(data.frame(Message = "No statistical test selected."))
    }
    
    res <- analysisResults()
    req(res)
    
    if (input$statTest %in% c("anova", "welch.anova", "kruskal")) {
      if (input$statTest == "anova") {
        result <- summary(res$statsResult)[[1]]
        if ("Pr(>F)" %in% colnames(result)) {
          data.frame(
            Term = rownames(result),
            `F value` = round(result$`F value`, 3),
            `Pr(>F)` = signif(result$`Pr(>F)`, 3),
            row.names = NULL
          )
        } else {
          data.frame(
            Statistic = names(res$statsResult$statistic),
            Value = round(as.numeric(res$statsResult$statistic), 3),
            `p-value` = signif(res$statsResult$p.value, 3),
            row.names = NULL
          )
        }
      } else if (input$statTest == "welch.anova") {
        result <- res$statsResult
        data.frame(
          `Statistic` = "F",
          `Value` = round(result$statistic, 3),
          `p-value` = signif(result$p.value, 3),
          `DF1` = round(result$parameter[1], 3),
          `DF2` = round(result$parameter[2], 3),
          row.names = NULL
        )
      } else if (input$statTest == "kruskal") {
        result <- res$statsResult
        data.frame(
          `Statistic` = "Chi-squared",
          `Value` = round(result$statistic, 3),
          `p-value` = signif(result$p.value, 3),
          `DF` = result$parameter,
          row.names = NULL
        )
      }
    } else {
      if (input$statTest == "wilcox") {
        result <- res$statsResult$p.value
        df <- as.data.frame(as.table(result))
        colnames(df) <- c("Group 1", "Group 2", "p-value")
        df <- df[!is.na(df$`p-value`), ]
        df$`p-value` <- signif(df$`p-value`, 3)
        return(df)
      } else {
        data.frame(
          Statistic = names(res$statsResult$statistic),
          Value = round(as.numeric(res$statsResult$statistic), 3),
          `p-value` = signif(res$statsResult$p.value, 3),
          conf.int.lower = round(res$statsResult$conf.int[1], 3),
          conf.int.upper = round(res$statsResult$conf.int[2], 3),
          row.names = NULL
        )
      }
    }
  })
  
  # Download Handlers
  output$downloadPlot <- downloadHandler(
    filename = function() { "growth_plot.png" },
    content = function(file) {
      res <- analysisResults()
      req(res)
      
      p <- ggplot(res$summarizedData, aes(x = time, y = mean_growth, color = group)) +
        geom_point(size = 3) +  
        geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, linewidth = 1) +  
        geom_errorbar(aes(ymin = mean_growth - stderr, ymax = mean_growth + stderr), width = 0.2) +
        labs(
          title = input$plotTitle,
          x = input$xLabel,
          y = input$yLabel,
          color = "Group"
        ) +
        theme_light()
      
      if (nrow(res$annotations) > 0) {
        p <- p + geom_text(
          data = res$annotations, 
          aes(x = x, y = y, label = label), 
          size = 6, 
          vjust = -0.5,
          inherit.aes = FALSE
        )
      }
      
      ggsave(filename = file, plot = p, dpi = 800, width = 12, height = 8)
    }
  )
  
  output$downloadResults <- downloadHandler(
    filename = function() { "analysis_results.txt" },
    content = function(file) {
      res <- analysisResults()
      req(res)
      
      sink(file)
      cat("Statistical Test Results:\n")
      if (input$statTest == "none") {
        cat("No statistical test selected.\n")
      } else if (input$statTest %in% c("anova", "welch.anova", "kruskal")) {
        if (input$statTest == "anova") {
          print(summary(res$statsResult))
        } else {
          print(res$statsResult)
        }
      } else if (input$statTest == "wilcox") {
        print(res$statsResult)
      } else {
        print(res$statsResult)
      }
      cat("\nPolynomial Regression Equations:\n")
      print(res$regressionEquations)
      cat("\nDoubling Time and Growth Rate (Exponential Model):\n")
      print(res$doublingGrowth)
      sink()
    }
  )
}

# Run the Shiny App
shinyApp(ui = ui, server = server)




