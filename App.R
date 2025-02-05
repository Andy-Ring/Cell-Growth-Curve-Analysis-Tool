# Load Required Libraries
library(shiny)
library(ggplot2)
library(dplyr)
library(car)  
library(stats)
library(ggpubr)  
library(multcomp)  # For multiple comparisons, if needed

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

# UI Definition
ui <- fluidPage(
  titlePanel("Cell Growth Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      downloadButton("downloadTemplate", "Download Template"),
      fileInput("data", "Upload Growth Data (CSV)", accept = ".csv"),
      actionButton("recommendTest", "Recommend Best Test"),
      textOutput("recommendedTest"),
      selectInput("statTest", "Select Statistical Test", 
                  choices = c(
                    "T-test" = "t.test", 
                    "ANOVA" = "anova", 
                    "Welch's ANOVA" = "welch.anova", 
                    "Kruskal-Wallis" = "kruskal",
                    "Pairwise Wilcoxon" = "wilcox"
                  )),
      actionButton("goButton", "Run Analysis"),
      downloadButton("downloadPlot", "Download High-Res Plot"),
      downloadButton("downloadResults", "Download Results")
    ),
    
    mainPanel(
      plotOutput("growthPlot"),
      tableOutput("statsTable"),
      verbatimTextOutput("regressionEquations")  
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
    req(growthData)  # Ensure data is available
    
    groups <- unique(growthData$group)
    recommendation <- ""
    
    if (length(groups) == 2) {
      # Check Normality
      normTest <- tryCatch(shapiro.test(growthData$growth)$p.value, error = function(e) NA)
      # Check Homogeneity of Variance
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
      # More than 2 groups
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
  
  # Run Analysis upon Clicking "Run Analysis"
  observeEvent(input$goButton, {
    growthData <- data()
    req(growthData)  # Ensure data is available
    
    # Convert group to factor (redundant if already done in data())
    growthData$group <- as.factor(growthData$group)
    
    # Summarize Data
    summarizedData <- growthData %>%
      group_by(time, group) %>%
      summarise(
        mean_growth = mean(growth, na.rm = TRUE),
        stderr = sd(growth, na.rm = TRUE) / sqrt(n()),
        .groups = 'drop'
      )
    
    # Generate Regression Equations
    regressionEquations <- summarizedData %>%
      group_by(group) %>%
      do({
        # Ensure enough data points for polynomial regression
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
    
    # Render Regression Equations
    output$regressionEquations <- renderPrint({
      cat("Regression Equations:\n")
      if(nrow(regressionEquations) == 0){
        cat("No regression equations available.\n")
      } else {
        apply(regressionEquations, 1, function(row) {
          cat(row["group"], ":", row["equation"], "\n")
        })
      }
    })
    
    # Perform Selected Statistical Test
    statsResult <- switch(input$statTest,
                          "t.test" = {
                            if(length(levels(growthData$group)) == 2){
                              t.test(growth ~ group, data = growthData)
                            } else {
                              showModal(modalDialog(
                                title = "Statistical Test Error",
                                "T-test is selected but there are more than 2 groups. Please select an appropriate test.",
                                easyClose = TRUE
                              ))
                              return(NULL)
                            }
                          },
                          "anova" = {
                            if(length(levels(growthData$group)) > 1){
                              aov(growth ~ group, data = growthData)
                            } else {
                              showModal(modalDialog(
                                title = "Statistical Test Error",
                                "ANOVA requires more than one group.",
                                easyClose = TRUE
                              ))
                              return(NULL)
                            }
                          },
                          "welch.anova" = {
                            if(length(levels(growthData$group)) > 1){
                              oneway.test(growth ~ group, data = growthData, var.equal = FALSE)
                            } else {
                              showModal(modalDialog(
                                title = "Statistical Test Error",
                                "Welch's ANOVA requires more than one group.",
                                easyClose = TRUE
                              ))
                              return(NULL)
                            }
                          },
                          "kruskal" = {
                            if(length(levels(growthData$group)) > 1){
                              kruskal.test(growth ~ group, data = growthData)
                            } else {
                              showModal(modalDialog(
                                title = "Statistical Test Error",
                                "Kruskal-Wallis Test requires more than one group.",
                                easyClose = TRUE
                              ))
                              return(NULL)
                            }
                          },
                          "wilcox" = {
                            if(length(levels(growthData$group)) > 1){
                              pairwise.wilcox.test(growthData$growth, growthData$group, p.adjust.method = "bonferroni")
                            } else {
                              showModal(modalDialog(
                                title = "Statistical Test Error",
                                "Pairwise Wilcoxon Test requires more than one group.",
                                easyClose = TRUE
                              ))
                              return(NULL)
                            }
                          }
    )
    
    # Check if statsResult is NULL due to earlier errors
    if(is.null(statsResult)){
      return(NULL)
    }
    
    # Render Statistical Test Results
    output$statsTable <- renderTable({
      if (input$statTest %in% c("anova", "welch.anova", "kruskal")) {
        if (input$statTest == "anova") {
          result <- summary(statsResult)[[1]]
          # Check if it's ANOVA or another type
          if ("Pr(>F)" %in% colnames(result)) {
            data.frame(
              Term = rownames(result),
              `F value` = round(result$`F value`, 3),
              `Pr(>F)` = signif(result$`Pr(>F)`, 3),
              row.names = NULL
            )
          } else {
            data.frame(
              Statistic = names(statsResult$statistic),
              Value = round(as.numeric(statsResult$statistic), 3),
              `p-value` = signif(statsResult$p.value, 3),
              row.names = NULL
            )
          }
        } else if (input$statTest == "welch.anova") {
          result <- statsResult
          data.frame(
            `Statistic` = "F",
            `Value` = round(result$statistic, 3),
            `p-value` = signif(result$p.value, 3),
            `DF1` = round(result$parameter[1], 3),
            `DF2` = round(result$parameter[2], 3),
            row.names = NULL
          )
        } else if (input$statTest == "kruskal") {
          result <- statsResult
          data.frame(
            `Statistic` = "Chi-squared",
            `Value` = round(result$statistic, 3),
            `p-value` = signif(result$p.value, 3),
            `DF` = result$parameter,
            row.names = NULL
          )
        }
      } else {
        # For t-test and Wilcoxon
        if (input$statTest == "wilcox") {
          # Pairwise Wilcoxon
          result <- statsResult$p.value
          # Melt the matrix for display
          df <- as.data.frame(as.table(result))
          colnames(df) <- c("Group 1", "Group 2", "p-value")
          df <- df[!is.na(df$`p-value`), ]
          df$`p-value` <- signif(df$`p-value`, 3)
          return(df)
        } else {
          # T-test
          data.frame(
            Statistic = names(statsResult$statistic),
            Value = round(as.numeric(statsResult$statistic), 3),
            `p-value` = signif(statsResult$p.value, 3),
            conf.int.lower = round(statsResult$conf.int[1], 3),
            conf.int.upper = round(statsResult$conf.int[2], 3),
            row.names = NULL
          )
        }
      }
    }, rownames = FALSE)
    
    # --- Revised Section: Handle Pairwise Comparisons for Annotation ---
    pairwise_comparisons <- NULL
    if (input$statTest == "anova") {
      # Perform ANOVA and Tukey's HSD
      anova_model <- aov(growth ~ group, data = growthData)
      tukey_results <- TukeyHSD(anova_model)$group
      
      # Create a matrix to store p-values
      groups <- levels(growthData$group)
      pairwise_comparisons <- matrix(NA, nrow = length(groups), ncol = length(groups), 
                                     dimnames = list(groups, groups))
      
      # Fill the matrix with adjusted p-values from Tukey's HSD
      for (comparison in rownames(tukey_results)) {
        groups_pair <- unlist(strsplit(comparison, "-"))
        if (length(groups_pair) == 2) {
          g1 <- groups_pair[1]
          g2 <- groups_pair[2]
          p_adj <- tukey_results[comparison, "p adj"]
          pairwise_comparisons[g1, g2] <- p_adj
          pairwise_comparisons[g2, g1] <- p_adj  # Ensure symmetry
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
      # Determine a vertical offset based on the error bar height (adjust multiplier as needed)
      offset <- max(summarizedData$stderr, na.rm = TRUE) * 2  
      annotation_counter <- 0
      
      for (combo in group_combos) {
        group1 <- combo[1]
        group2 <- combo[2]
        # Check if the matrix contains a p-value for this pair
        if(all(c(group1, group2) %in% rownames(pairwise_comparisons))){
          p_val <- pairwise_comparisons[group1, group2]
        } else {
          p_val <- NA
        }
        if (!is.na(p_val)) {
          annotation_label <- significance_labels(p_val)
          annotation_counter <- annotation_counter + 1
          # Place each annotation at the middle of the time range; adjust vertical position so they stack
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
    
    # Render Growth Plot with Static Significance Annotations
    output$growthPlot <- renderPlot({
      req(summarizedData)  # Ensure data is available
      
      p <- ggplot(summarizedData, aes(x = time, y = mean_growth, color = group)) +
        geom_point(size = 3) +  
        geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, linewidth = 1) +  
        geom_errorbar(aes(ymin = mean_growth - stderr, ymax = mean_growth + stderr), width = 0.2) +  
        labs(x = "Time", y = "Mean Growth", color = "Group") +
        theme_light()
      
      # Add the significance annotations if any exist
      if (nrow(annotations) > 0) {
        p <- p + geom_text(
          data = annotations, 
          aes(x = x, y = y, label = label), 
          size = 6, 
          vjust = -0.5,
          inherit.aes = FALSE
        )
      }
      
      p
    })
    
    # Enable Downloading of Plot
    output$downloadPlot <- downloadHandler(
      filename = function() { "growth_plot.png" },
      content = function(file) {
        png(file, width = 1600, height = 1200, res = 150)
        p <- ggplot(summarizedData, aes(x = time, y = mean_growth, color = group)) +
          geom_point(size = 3) +  
          geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE, linewidth = 1) +  
          geom_errorbar(aes(ymin = mean_growth - stderr, ymax = mean_growth + stderr), width = 0.2) +  
          labs(x = "Time", y = "Mean Growth", color = "Group") +
          theme_light()
        
        if (nrow(annotations) > 0) {
          p <- p + geom_text(
            data = annotations, 
            aes(x = x, y = y, label = label), 
            size = 6, 
            vjust = -0.5,
            inherit.aes = FALSE
          )
        }
        print(p)
        dev.off()
      }
    )
    
    # Prepare and Enable Downloading of Results
    output$downloadResults <- downloadHandler(
      filename = function() { "analysis_results.txt" },
      content = function(file) {
        sink(file)
        cat("Statistical Test Results:\n")
        if (input$statTest %in% c("anova", "welch.anova", "kruskal")) {
          if (input$statTest == "anova") {
            print(summary(statsResult))
          } else {
            print(statsResult)
          }
        } else if (input$statTest == "wilcox") {
          print(statsResult)
        } else {
          print(statsResult)
        }
        cat("\nRegression Equations:\n")
        print(regressionEquations)
        sink()
      }
    )
    
  })
}

# Run the Shiny App
shinyApp(ui = ui, server = server)


