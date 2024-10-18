deconvolution <- function(days, real_deaths, num_reps = 100, bootstrap = FALSE, initial_t0 = NULL) {
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
    inf_to_death <- sample(1:max_duration, total_deaths, replace = TRUE, prob = probability_vector)
    initial_t0 <- death_days - inf_to_death
  }

  incidence_matrix <- matrix(0, nrow = max_days, ncol = num_reps)
  p_statistics <- numeric(num_reps)
  inf_to_death <- sample(1:max_duration, total_deaths, replace = TRUE, prob = probability_vector)

  lower_bounds <- matrix(0, nrow = max_days, ncol = num_reps)
  upper_bounds <- matrix(0, nrow = max_days, ncol = num_reps)

  no_improvement_count <- 0
  bootstrap_results <- matrix(0, nrow = max_days, ncol = num_reps)

  for (rep in seq_len(num_reps)) {
    simulated_deaths <- initial_t0 + inf_to_death
    daily_deaths <- tabulate(simulated_deaths, nbins = max_days)

    p_statistic <- sum((real_deaths - daily_deaths[1:length(real_deaths)]) ^ 2 / pmax(1, daily_deaths[1:length(real_deaths)]))
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

      new_p_statistic <- sum((real_deaths - daily_deaths[1:length(real_deaths)]) ^ 2 / pmax(1, daily_deaths[1:length(real_deaths)]));

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

    if (no_improvement_count >= patience_threshold) {
      message("Early stopping at iteration ", rep, " due to no improvement.")
      break
    }
  }

  result <- list(PValues = p_statistics, IncidenceTrajectory = incidence_matrix, InitialT0 = initial_t0)
  return(result)
}

# setwd(file.path(current_wd, "data"))                 # (Uncomment this to run)
# data <- read.table("engcov.txt", header = TRUE)      # (Uncomment this to run)
days <- data$julian[1:150]
real_deaths <- data$nhs[1:150]

result <- deconvolution(days, real_deaths, num_reps = 100, bootstrap = TRUE)
