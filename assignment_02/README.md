# Inferring Fatal Incidence Rates from COVID-19 Death Data

## Introduction

During the COVID-19 pandemic, understanding the dynamics of new infections was critical for managing the epidemic. The incidence of new infections indicates whether the epidemic is under control, and it is a more immediate measure than the number of confirmed cases, which often lag behind the actual infections. However, new infections are challenging to observe directly; instead, we can infer them from the distribution of COVID-19 deaths over time. 

This article outlines a method for estimating the fatal incidence rates from daily death data. By leveraging the distribution of time from infection to death, we can better understand how the epidemic evolved over time and assess the effectiveness of control measures.

## Concept Overview

The core idea behind this method involves modeling the relationship between infection times and death times using a statistical distribution, specifically a log-normal distribution. We assign each death a guessed time of infection and then simulate the corresponding death times by adding random draws from the infection-to-death distribution.

### Convolution and Deconvolution

To understand how we can derive the incidence from death data, we first need to discuss the concepts of convolution and deconvolution:

- **Convolution** is a mathematical operation that combines two functions to produce a third function, representing the amount of overlap between the two functions. In our context, if we denote the infection distribution as \( f(t) \) and the time from infection to death as \( g(t) \), the convolution \( (f * g)(t) \) gives us the distribution of deaths over time. 

  \[
  (f * g)(t) = \int_{-\infty}^{\infty} f(\tau) g(t - \tau) d\tau
  \]

  This operation allows us to model how infections at different times contribute to the number of deaths observed at a given time.

- **Deconvolution**, on the other hand, is the process of extracting one function from the convolution of two functions. In our case, we want to reverse the convolution process to infer the number of new infections based on the observed deaths. This involves estimating the original infection distribution \( f(t) \) given the death data and the known time-to-death distribution \( g(t) \).

## Simulation Methodology

### Step 1: Infection-to-Death Distribution

First, we need to determine the probabilities of infection-to-death durations. We utilize a log-normal distribution characterized by parameters derived from the ISARIC study, specifically \( \text{meanlog} = 3.152 \) and \( \text{sdlog} = 0.451 \). The log-normal distribution is appropriate for modeling positive skewed data, which aligns with our understanding of how long it takes from infection to death.

The log-normal probability density function is given by:

\[
g(t) = \frac{1}{t \sigma \sqrt{2\pi}} e^{-\frac{(\ln t - \mu)^2}{2\sigma^2}}
\]

where \( \mu = \text{meanlog} \) and \( \sigma = \text{sdlog} \).

### Step 2: Initial Guess for Infection Times

With the infection-to-death distribution defined, we can begin the simulation process:

1. **Generate Initial Infection Times**: For each of the fatalities recorded in the dataset, we will randomly sample from the infection-to-death distribution to create an initial guess for the infection times. This will yield a vector \( t_0 \) representing the guessed infection days corresponding to each death day.

### Step 3: Iterative Refinement

To improve our estimates, we will iteratively refine the infection times:

1. **Simulate Death Days**: For each proposed infection time in \( t_0 \), draw new values from the infection-to-death distribution to generate simulated death days. This process involves updating our estimates based on the current state of the infection times.

2. **Goodness of Fit**: We evaluate how well our simulated death counts match the observed counts using a modified Pearson statistic:

   \[
   P = \sum_{i} \frac{(d_i - d_{si})^2}{\max(1, d_{si})}
   \]

   where \( d_i \) is the actual number of deaths on day \( i \), and \( d_{si} \) is the number of deaths according to our simulation.

3. **Propose Adjustments**: For each infection time, we randomly propose moving it by a few days (either forward or backward). If the proposed change improves the fit \( P \), we accept it; otherwise, we maintain the original time. This random walk approach continues iteratively, refining our estimates.

## Conclusion

This simulation method effectively uses the statistical properties of the infection-to-death distribution to infer fatal incidence rates from COVID-19 death data. By understanding the convolution of distributions and applying iterative refinement techniques, we can derive valuable insights into the epidemic's progression.

As we observe changes in the death data over time, we gain a clearer understanding of the underlying infection dynamics, which ultimately aids public health officials in making informed decisions about control measures and resource allocation.

In summary, by employing statistical modeling and simulation, we can bridge the gap between observed deaths and unobserved infections, thus gaining a more comprehensive view of the COVID-19 pandemic's trajectory.

## References

1. ISARIC Study
2. Log-Normal Distribution