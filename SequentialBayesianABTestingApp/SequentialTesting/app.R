# Load libraries required for the Shiny app, Bayesian analysis, and plotting
library(bayesAB)
library(shiny)
library(ggplot2)

# Custom function: Updates the posterior distribution using the prior distribution and observed data
update_normal <- function(prior, data) {
  # Calculate new mean and variance for the updated posterior distribution
  new_mean <- (prior$mean / prior$var + sum(data) / (length(data) * prior$var)) / (1 / prior$var + length(data) / (length(data) * prior$var))
  new_var <- 1 / (1 / prior$var + length(data) / (length(data) * prior$var))
  
  # Compute the likelihood values for each data point
  likelihood <- dnorm(data, mean = mean(data), sd = sqrt(prior$var))
  print(paste("Likelihood values:", paste(likelihood, collapse = ", ")))
  
  # Return the updated posterior distribution as a list
  return(list(mean = new_mean, var = new_var))
}

# Custom function: Calculates the probability that the mean of version A is greater than the mean of version B
prob_superiority <- function(posterior_A, posterior_B) {
  return(pnorm((posterior_A$mean - posterior_B$mean) / sqrt(posterior_A$var + posterior_B$var)))
}

# Create the Shiny app user interface layout with inputs and outputs
ui <- fluidPage(
  titlePanel("Bayesian AB Testing with Sequential Testing"),
  sidebarLayout(
    sidebarPanel(
      numericInput("n", "Maximum Observations:", value = 1000, min = 1, max = 10000),
      sliderInput("prior_mean", "Initial Average Effect:", min = -1, max = 1, value = 0.2, step = 0.01),
      sliderInput("prior_sd", "Initial Effect Variation:", min = 0.01, max = 1, value = 0.1, step = 0.01),
      numericInput("min_effect_size", "Minimum Acceptable Effect Size:", value = 0.2, min = 0, max = 1, step = 0.01),
      numericInput("stop_rule", "Confidence Level (as a percentage, i.e., 99%):", value = 0.95, min = 0.5, max = 1, step = 0.01),
      numericInput("mean_A", "Average Performance of Version A:", value = 0.4, min = -1, max = 1, step = 0.01),
      numericInput("sd_A", "Performance Variability of Version A:", value = 0.15, min = 0.01, max = 1, step = 0.01),
      numericInput("mean_B", "Average Performance of Version B:", value = 0.5, min = -1, max = 1, step = 0.01),
      numericInput("sd_B", "Performance Variability of Version B:", value = 0.15, min = 0.01, max = 1, step = 0.01),
      actionButton("goButton", "Start Test")
    ),
    mainPanel(
      plotOutput("plot"),
      plotOutput("thresholdPlot"),
      textOutput("stoppingInfo"),
      textOutput("stoppingStepInfo"),
      textOutput("summaryInfo"),
      textOutput("probSuperiorityInfo")
    )
  )
)

# Define the Shiny app server function, which runs the Bayesian AB testing and generates the output plots and summaries
server <- function(input, output) {
  observeEvent(input$goButton, {
    
    # Set the random seed for reproducibility of the results
    set.seed(123)
    
    # Initialize the prior distribution using the user inputs for mean and variance
    prior <- list(mean = input$prior_mean, var = input$prior_sd^2)
    
    # Initialize the posterior distribution for versions A and B using the prior distribution
    posterior_A <- prior
    posterior_B <- prior
    
    # Initialize the sample sizes for versions A and B
    n_A <- 0
    n_B <- 0
    
    # Initialize the observed data for versions A and B
    data_A <- numeric()
    data_B <- numeric()
    
    # Initialize a data frame to store the stopping threshold values at each step
    threshold_values <- data.frame(step = integer(), value = numeric())
    
    # Perform Bayesian sequential testing by collecting data until the stopping rule is met
    while (TRUE) {
      
      # Choose which version (A or B) to sample from based on the current step
      if ((n_A + n_B) %% 2 == 0) {
        n_A <- n_A + 1
        data_A[n_A] <- rnorm(1, mean = input$mean_A, sd = input$sd_A)
        posterior_A <- update_normal(posterior_A, data_A)
      } else {
        n_B <- n_B + 1
        data_B[n_B] <- rnorm(1, mean = input$mean_B, sd = input$sd_B)
        posterior_B <- update_normal(posterior_B, data_B)
      }
      
      # Calculate the probability that version A is better than version B using the prob_superiority custom function
      prob_A_better <- prob_superiority(posterior_A, posterior_B)
      
      # Calculate the difference in means between the two treatments based on their posterior distributions
      mean_diff <- posterior_A$mean - posterior_B$mean
      
      # Add the current probability value to the threshold_values data frame for plotting
      threshold_values <- rbind(threshold_values, data.frame(step = n_A + n_B, value = prob_A_better))
      
      # Check if the stopping rule is met based on the probability of superiority and the user-defined minimum acceptable effect size
      if (mean_diff >= input$min_effect_size) {
        if (prob_A_better > input$stop_rule | prob_A_better < 1 - input$stop_rule) {
          break
        }
      } else {
        if (prob_A_better > input$stop_rule & mean_diff >= input$min_effect_size | 
            prob_A_better < 1 - input$stop_rule & mean_diff <= -input$min_effect_size) {
          break
        }
      }
      
      # Print statements for debugging purposes, displaying the current state of the analysis
      print(paste0("n_A = ", n_A))
      print(paste0("data_A = ", data_A))
      print(paste0("posterior_A = ", posterior_A))
      
    }
    
    # Display stopping step and threshold in the Shiny app output
    output$stoppingStepInfo <- renderText({
      paste0("Stopping Step: ", n_A + n_B, "\n",
             "Threshold: ", input$stop_rule)
    })
    
    # Display sample sizes, posterior means and variances, and probability of superiority in the Shiny app output
    output$summaryInfo <- renderText({
      paste0("Sample sizes:\n",
             "Version A: ", n_A, "\n",
             "Version B: ", n_B, "\n",
             "Posterior means:\n",
             "Version A: ", round(posterior_A$mean, 4), "\n",
             "Version B: ", round(posterior_B$mean, 4), "\n",
             "Posterior variances:\n",
             "Version A: ", round(posterior_A$var, 4), "\n",
             "Version B: ", round(posterior_B$var, 4))
    })
    
    # Display probability of superiority in the Shiny app output
    output$probSuperiorityInfo <- renderText({
      paste0("Probability of Superiority:\n",
             "P(A > B) = ", round(prob_A_better, 4))
    })
    
    # Display minimum acceptable effect size in the Shiny app output
    output$minEffectSizeInfo <- renderText({
      paste0("Minimum acceptable effect size: ", input$min_effect_size)
    })
    
    # Generate the posterior distribution plot comparing versions A and B
    output$plot <- renderPlot({
      data_frame <- data.frame(version = c(rep("A", length(data_A)), rep("B", length(data_B))),
                               value = c(data_A, data_B))
      p <- ggplot(data_frame, aes(x = value, fill = version)) +
        geom_density(alpha = 0.4) +
        labs(title = "Posterior Distribution: Updated Beliefs About Model Performance for Versions A and B",
             x = "Effect Difference",
             y = "Probability Density") +
        scale_x_continuous(limits = c(min(c(input$mean_A - 3 * input$sd_A, input$mean_B - 3 * input$sd_B)), max(c(input$mean_A + 3 * input$sd_A, input$mean_B + 3 * input$sd_B)))) +
        scale_y_continuous(limits = c(0, max(density(c(data_A, data_B))$y, na.rm = TRUE, finite = TRUE) * 1.5)) +
        theme_minimal()
      
      # Add label for the probability of superiority
      p <- p + annotate("text", x = max(data_frame$value), y = max(density(c(data_A, data_B))$y, na.rm = TRUE, finite = TRUE),
                        label = paste("P(A > B) =", round(prob_A_better, 4)),
                        hjust = 1, vjust = 1, size = 4)
      p
    })
    
    # Generate the stopping threshold plot with the vertical line indicating
    # when the stopping rule was met
    output$thresholdPlot <- renderPlot({
      # Create a ggplot object to plot the threshold values over time
      ggplot(threshold_values, aes(x = step, y = value)) +
        # Add a line plot of the threshold values
        geom_line() +
        # Add red dashed lines for the stopping rule upper and lower thresholds
        geom_hline(yintercept = input$stop_rule, linetype = "dashed", color = "red") +
        geom_hline(yintercept = 1 - input$stop_rule, linetype = "dashed", color = "red") +
        # Add a blue dashed vertical line at the step where the stopping rule was met
        geom_vline(xintercept = n_A + n_B, linetype = "dashed", color = "blue") + 
        # Add plot labels and use a minimal theme
        labs(title = "Stopping Threshold Over Time: Monitoring Decision Confidence",
             x = "Steps",
             y = "Probability of Superiority") +
        theme_minimal()
    })
  })
}

# Run the Shiny app with the defined user interface and server function
shinyApp(ui = ui, server = server)
