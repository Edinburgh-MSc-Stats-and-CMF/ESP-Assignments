# -------------------------------- TEAM CONTRIBUTIONS ------------------------ #
# - Sagar Udasi: S2606876
#            40% - Developed the main deconvolution algorithm and
#                        optimized the bootstrap functionality.
# - Xu Guo: S2743905
#            30% - Created visualization functions for plotting real vs.
#                        simulated deaths and P statistics.
# - Xie Ying: S2510813
#            30% - Assisted in topic research and understanding,
#                        annotations, and contributed to the overall testing
#                        and validation of the code.

# -------------------------- PROBLEM STATEMENT AND APPROACH ------------------ #
# This code models the distribution of COVID-19 deaths over time using a
# deconvolution approach. The goal is to estimate the number of infections that
# resulted in reported deaths, accounting for the time lag between infection
# and death. The code simulates death trajectories based on historical data and
# generates visualizations to assess the fit of simulated data against
# real-world observations.

# --------------------------------- CODE ------------------------------------- #
library(ggplot2)

if (!dir.exists("plots")) dir.create("plots", showWarnings = FALSE)

# This function generates and saves a plot comparing real deaths and simulated
# deaths with confidence intervals.
# Inputs: days (numeric), real_deaths (numeric), simulated_deaths (numeric),
#         lower_bounds (numeric), upper_bounds (numeric), iteration (integer)
plot_real_vs_simulated <- function(days,
                                   real_deaths,
                                   simulated_deaths,
                                   lower_bounds,
                                   upper_bounds,
                                   iteration) {
  plot_data <- data.frame(
    Day = days,
    RealDeaths = real_deaths,
    SimulatedDeaths = simulated_deaths[1:length(real_deaths)],
    LowerBound = lower_bounds,
    UpperBound = upper_bounds
  )

  plot <- ggplot(plot_data, aes(x = Day)) +
    geom_line(aes(y = RealDeaths, color = "Real Deaths"), linewidth = 1) +
    geom_line(aes(y = SimulatedDeaths, color = "Simulated Deaths"),
              linewidth = 1) +
    geom_ribbon(aes(ymin = LowerBound, ymax = UpperBound),
                fill = "blue",
                alpha = 0.2) +
    labs(title = paste("Iteration:", iteration),
         x = "Day",
         y = "Deaths") +
    scale_color_manual(values = c(
      "Real Deaths" = "red",
      "Simulated Deaths" = "blue"
    ))

  plot <- plot + geom_vline(xintercept = 84,
                            linetype = "dashed",
                            color = "black") +
    annotate(
      "text",
      x = 84,
      y = max(real_deaths, na.rm = TRUE),
      label = "UK Lockdown\nMarch 24, 2020",
      hjust = -0.1,
      color = "black"
    )

  ggsave(
    file.path("plots", paste(
      "iteration_plot_", iteration, ".png", sep = ""
    )),
    plot = plot,
    width = 10,
    height = 6,
    bg = "white"
  )
}

# Function generates and saves a plot for P statistic values over iterations.
# Inputs: p_values (numeric vector)
plot_p_statistic <- function(p_values) {
  p_statistic_data <- data.frame(Iteration = seq_along(p_values),
                                 PStatistic = p_values)

  p_statistic_plot <- ggplot(p_statistic_data, aes(x = Iteration,
                                                   y = PStatistic)) +
    geom_line(color = "blue") +
    geom_point(color = "blue") +
    labs(x = "Iteration", y = "P Statistic", 
         title = "P Statistic Over Iterations")

  ggsave(
    file.path("plots", "p_statistic_plot.png"),
    plot = p_statistic_plot,
    width = 10,
    height = 6,
    bg = "white"
  )
}

# MAIN DECONVOLUTION FUNCTION:
# This function estimates the distribution of deaths over time using a
# deconvolution approach.
# Inputs: days (numeric vector), real_deaths (numeric vector),
#         num_reps (integer), bootstrap (boolean),
#         initial_t0 (optional numeric vector)
# Outputs: list containing PValues, IncidenceTrajectory, and InitialT0
deconvolution <- function(days,
                          real_deaths,
                          num_reps = 100,
                          bootstrap = FALSE,
                          initial_t0 = NULL) {
  log_mean <- 3.152
  log_sd <- 0.451
  max_duration <- 80
  max_days <- 150
  tolerance <- 1e-5
  patience_threshold <- 20

  probability_vector <- dlnorm(1:max_duration, log_mean, log_sd)
  probability_vector <- probability_vector / sum(probability_vector)

  total_deaths <- sum(real_deaths)

  if (is.null(initial_t0)) {
    death_days <- rep(days, real_deaths)
    inf_to_death <- sample(1:max_duration,
                           total_deaths,
                           replace = TRUE,
                           prob = probability_vector)
    initial_t0 <- death_days - inf_to_death
  }

  incidence_matrix <- matrix(0, nrow = max_days, ncol = num_reps)
  p_statistics <- numeric(num_reps)
  inf_to_death <- sample(1:max_duration,
                         total_deaths,
                         replace = TRUE,
                         prob = probability_vector)

  lower_bounds <- matrix(0, nrow = max_days, ncol = num_reps)
  upper_bounds <- matrix(0, nrow = max_days, ncol = num_reps)

  no_improvement_count <- 0
  bootstrap_results <- matrix(0, nrow = max_days, ncol = num_reps)

  for (rep in seq_len(num_reps)) {
    simulated_deaths <- initial_t0 + inf_to_death
    daily_deaths <- tabulate(simulated_deaths, nbins = max_days)

    p_statistic <- sum((real_deaths - daily_deaths[1:length(real_deaths)]) ^
                         2 / pmax(1, daily_deaths[1:length(real_deaths)]))
    p_statistics[rep] <- p_statistic

    incidence_matrix[, rep] <- tabulate(initial_t0, nbins = max_days)

    lower_bounds[, rep] <- pmax(0, daily_deaths - 1.96 * sqrt(daily_deaths))
    upper_bounds[, rep] <- daily_deaths + 1.96 * sqrt(daily_deaths)

    indices <- sample(total_deaths)
    steps <- sample(c(-8, -4, -2, -1, 1, 2, 4, 8), total_deaths, replace = TRUE)

    for (i in indices) {
      current_t0 <- initial_t0[i]
      step <- steps[i]
      new_t0 <- current_t0 + step

      if (new_t0 < 1 || new_t0 > max_days)
        next

      old_death_day <- simulated_deaths[i]
      new_death_day <- new_t0 + inf_to_death[i]

      if (old_death_day != new_death_day) {
        daily_deaths[old_death_day] <- daily_deaths[old_death_day] - 1
        daily_deaths[new_death_day] <- daily_deaths[new_death_day] + 1
      }

      new_p_statistic <- sum((real_deaths - daily_deaths[1:length(real_deaths)])
                             ^ 2 / pmax(1, daily_deaths[1:length(real_deaths)]))

      if (new_p_statistic < p_statistic) {
        initial_t0[i] <- new_t0
        simulated_deaths[i] <- new_death_day
        p_statistic <- new_p_statistic
        no_improvement_count <- 0
      } else {
        daily_deaths[old_death_day] <- daily_deaths[old_death_day] + 1
        daily_deaths[new_death_day] <- daily_deaths[new_death_day] - 1
        no_improvement_count <- no_improvement_count + 1
      }
    }

    if (bootstrap) {
      bootstrap_results[, rep] <- daily_deaths[1:max_days]
    }

    if (rep == 1 || rep %% 5 == 0) {
      plot_real_vs_simulated(days,
                             real_deaths,
                             daily_deaths,
                             lower_bounds[, rep],
                             upper_bounds[, rep],
                             rep)
    }

    if (no_improvement_count >= patience_threshold) {
      message("Early stopping at iteration ", rep, " due to no improvement.")
      plot_real_vs_simulated(days,
                             real_deaths,
                             daily_deaths,
                             lower_bounds[, rep],
                             upper_bounds[, rep],
                             rep)
      break
    }
  }

  plot_p_statistic(p_statistics)

  result <- list(
    PValues = p_statistics,
    IncidenceTrajectory = incidence_matrix,
    InitialT0 = initial_t0
  )
  return(result)
}

# MAIN DRIVER CODE
# Load the dataset containing COVID-19 death counts and corresponding days.

# setwd(file.path(current_wd(), "data"))             # (Uncomment this to run)
data <- read.table("engcov.txt", header = TRUE)      # (Uncomment this to run)
days <- data$julian[1:150]
real_deaths <- data$nhs[1:150]

result <- deconvolution(days, real_deaths, num_reps = 100, bootstrap = TRUE)
