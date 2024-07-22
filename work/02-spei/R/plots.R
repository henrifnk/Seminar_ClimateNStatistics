library(tidyverse)
library(ggpubr)
library(zoo)
library(latex2exp)

climate_data <-  readRDS("work/02-spei/data/climate_data.rds")

time_scales <- c(1, 3, 6, 12, 18, 24)

theme_set(theme_bw())

# Function to estimate the parameters of the 3-parameter log-logistic distribution using L-Moments
L_moments <- function(data, col) {
  data <- data[!is.na(data[[col]]), ]
  
  N <- nrow(data)
  i <- 1:N
  freq_est <- (i - 0.35) / N
  w_0 <- 1 / N * sum((1 - freq_est) ^ 0 * sort(data[[col]]))
  w_1 <- 1 / N * sum((1 - freq_est) ^ 1 * sort(data[[col]]))
  w_2 <- 1 / N * sum((1 - freq_est) ^ 2 * sort(data[[col]]))
  w_3 <- 1 / N * sum((1 - freq_est) ^ 3 * sort(data[[col]]))
  beta_est <- (2 * w_1 - w_0) / (6 * w_1 - w_0 - 6 * w_2)
  alpha_est <- ((w_0 - 2 * w_1) * beta_est) / (gamma(1 + 1 / beta_est) * gamma(1 - 1 / beta_est))
  gamma_est <- w_0 - alpha_est * gamma(1 + 1 / beta_est) * gamma(1 - 1 / beta_est)
  return(list(alpha = alpha_est, beta = beta_est, gamma = gamma_est))
}

# density function of the 3-parameter log-logistic distribution
f <- function(x, alpha, beta, gamma) {
  beta / alpha * ((x - gamma) / alpha) ^ (beta - 1) * (1 + ((x - gamma) / alpha) ^ beta) ^ (-2)
}

# distribution function of the 3-parameter log-logistic distribution
F_loglogistic <- function(x, alpha, beta, gamma) {
  if (alpha > gamma) {
    ifelse(x <= gamma, 0, 1 / (1 + (alpha / (x - gamma)) ^ beta))
  } else {
    ifelse(x >= gamma, 1, 1 / (1 + (alpha / (x - gamma)) ^ beta))
  }
}

reference_period <- climate_data %>%
  filter(Jahr >= 1981 & Jahr <= 2010)

## Model fitting Plot
i <- 0
plots <- list()

for (time_scale in time_scales) {
  i <- i + 1
  
  L_moments_values <- L_moments(reference_period, paste0("D_accum_", time_scale))
  
  p <- reference_period %>%
    ggplot(aes(x = !!sym(paste0("D_accum_", time_scale)))) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 30,
                   fill = "grey80",
                   color = "black",
                   alpha = 0.7) +
    stat_function(fun = f,
                  aes(color = "Log-Logistic"),
                  linewidth = 1,
                  args = L_moments_values) +
    stat_function(fun = dnorm,
                  aes(color = "Normal"),
                  linewidth = 1,
                  args = list(
                    mean = mean(reference_period[[paste0("D_accum_", time_scale)]], na.rm = TRUE),
                    sd = sd(reference_period[[paste0("D_accum_", time_scale)]], na.rm = TRUE))) +
    labs(x = "D", y = "Density") +
    ggtitle(ifelse(time_scale == 1, "1 Month", paste(time_scale, "Months"))) +
    scale_color_manual(values = c("Normal" = "royalblue",
                                  "Log-Logistic" = "orange")) +
    theme(legend.title = element_blank())
  
  plots[[i]] <- p
}

modelfitting_plot <- ggarrange(plotlist = plots,
                               ncol = 3,
                               nrow = 2,
                               common.legend = TRUE,
                               legend = "bottom")

## Wasserstein distance

wasserstein_distances <- matrix(NA, nrow = length(time_scales), ncol = 3)
colnames(wasserstein_distances) <- c("Time Scale", "Normal", "Log-Logistic")

i <- 0
for (time_scale in time_scales) {
  i <- i + 1
  
  L_moments_values <- L_moments(reference_period, paste0("D_accum_", time_scale))
  
  empirical_cdf <- ecdf(reference_period[[paste0("D_accum_", time_scale)]])
  theoretical_cdf_normal <- function(x)
    pnorm(
      x,
      mean = mean(reference_period[[paste0("D_accum_", time_scale)]], na.rm = TRUE),
      sd = sd(reference_period[[paste0("D_accum_", time_scale)]], na.rm = TRUE)
    )
  theoretical_cdf_loglogistic <- function(x)
    F_loglogistic(x,
                  L_moments_values$alpha,
                  L_moments_values$beta,
                  L_moments_values$gamma)
  
  abs_diff_cdf_normal <- function(x)
    abs(empirical_cdf(x) - theoretical_cdf_normal(x))
  abs_diff_cdf_loglogistic <- function(x)
    abs(empirical_cdf(x) - theoretical_cdf_loglogistic(x))
  
  wasserstein_distance_normal <- integrate(
    abs_diff_cdf_normal,
    lower = -Inf,
    upper = Inf,
    subdivisions = 10000
  )$value
  wasserstein_distance_loglogistic <- integrate(
    abs_diff_cdf_loglogistic,
    lower = -Inf,
    upper = Inf,
    subdivisions = 10000
  )$value
  
  wasserstein_distances[i, ] <- c(
    time_scale,
    round(wasserstein_distance_normal, 4),
    round(wasserstein_distance_loglogistic, 4)
  )
}

## SPEI Bavaria Plot

long_data <- climate_data %>%
  dplyr::select(Jahr, Monat, Datum, spei1, spei3, spei6, spei12, spei18, spei24) %>%
  pivot_longer(cols = starts_with("spei"),
               names_to = "time_scale",
               values_to = "spei") %>%
  mutate(time_scale = factor(
    time_scale,
    levels = c("spei1", "spei3", "spei6", "spei12", "spei18", "spei24"),
    labels = c(
      TeX("$SPEI^{(1)}$"),
      TeX("$SPEI^{(3)}$"),
      TeX("$SPEI^{(6)}$"),
      TeX("$SPEI^{(12)}$"),
      TeX("$SPEI^{(18)}$"),
      TeX("$SPEI^{(24)}$")
    )
  )) %>%
  filter(!is.na(spei))

spei_plot <- long_data %>%
  filter(Jahr > 1920) %>%
  ggplot(aes(x = Datum, y = spei)) +
  geom_area() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year", y = "SPEI") +
  facet_wrap(
    ~ time_scale,
    ncol = 1,
    scales = "free_y",
    strip.position = "left",
    labeller = label_parsed
  ) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    panel.spacing = unit(0.1, "lines"),
    axis.title.y = element_blank()
  )

## SPEI12 Extreme Droughts Plot

extreme_droughts <- climate_data %>%
  filter(Jahr > 1920) %>%
  mutate(extreme = ifelse(spei12 < -2, 1, 0)) %>%
  # for 10 years sum the number of extreme droughts
  mutate(extreme_sum = rollapply(extreme, 120, sum, fill = NA, align = "right") / 120) %>%
  filter(!is.na(extreme_sum))


extreme_droughts_plot <- ggplot(extreme_droughts, aes(x = Datum, y = extreme_sum)) +
  geom_line() +
  labs(x = "Year", y = TeX("$prop_{i,120}^{(12)}$")) +
  coord_fixed(ratio = 40000)

## Extreme Droughts for different time scales Plot

extreme_droughts_time_scales_plot <-
  climate_data %>%
  filter(Jahr > 1920) %>%
  mutate(
    extreme1 = ifelse(spei1 < -2, 1, 0),
    extreme3 = ifelse(spei3 < -2, 1, 0),
    extreme6 = ifelse(spei6 < -2, 1, 0),
    extreme12 = ifelse(spei12 < -2, 1, 0),
    extreme18 = ifelse(spei18 < -2, 1, 0),
    extreme24 = ifelse(spei24 < -2, 1, 0)
  ) %>%
  mutate(
    extreme_sum1 = rollapply(extreme1, 120, sum, fill = NA, align = "right") / 120,
    extreme_sum3 = rollapply(extreme3, 120, sum, fill = NA, align = "right") / 120,
    extreme_sum6 = rollapply(extreme6, 120, sum, fill = NA, align = "right") / 120,
    extreme_sum12 = rollapply(extreme12, 120, sum, fill = NA, align = "right") / 120,
    extreme_sum18 = rollapply(extreme18, 120, sum, fill = NA, align = "right") / 120,
    extreme_sum24 = rollapply(extreme24, 120, sum, fill = NA, align = "right") / 120
  ) %>%
  pivot_longer(
    cols = starts_with("extreme_sum"),
    names_to = "time_scale",
    values_to = "extreme_sum"
  ) %>%
  mutate(time_scale = factor(
    time_scale,
    levels = c(
      "extreme_sum1",
      "extreme_sum3",
      "extreme_sum6",
      "extreme_sum12",
      "extreme_sum18",
      "extreme_sum24"
    ),
    labels =
      c(
        TeX("$prop_{i,120}^{(1)}$"),
        TeX("$prop_{i,120}^{(3)}$"),
        TeX("$prop_{i,120}^{(6)}$"),
        TeX("$prop_{i,120}^{(12)}$"),
        TeX("$prop_{i,120}^{(18)}$"),
        TeX("$prop_{i,120}^{(24)}")
      )
  )) %>%
  filter(!is.na(extreme_sum)) %>%
  ggplot(aes(x = Datum, y = extreme_sum)) +
  geom_line() +
  labs(x = "Year") +
  facet_wrap(
    ~ time_scale,
    ncol = 1,
    strip.position = "left",
    labeller = label_parsed
  ) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    panel.spacing = unit(0.1, "lines"),
    axis.title.y = element_blank()
  )

## Extreme Droughts for different periods Plot

extreme_droughts_periods_plot <- climate_data %>%
  filter(Jahr > 1920) %>%
  mutate(extreme = ifelse(spei12 < -2, 1, 0)) %>%
  mutate(
    extreme_sum1 = rollapply(extreme, 12, sum, fill = NA, align = "right") / 12,
    extreme_sum5 = rollapply(extreme, 60, sum, fill = NA, align = "right") / 60,
    extreme_sum10 = rollapply(extreme, 120, sum, fill = NA, align = "right") / 120,
    extreme_sum15 = rollapply(extreme, 180, sum, fill = NA, align = "right") / 180,
    extreme_sum20 = rollapply(extreme, 240, sum, fill = NA, align = "right") / 240
  ) %>%
  pivot_longer(
    cols = starts_with("extreme_sum"),
    names_to = "period",
    values_to = "extreme_sum"
  ) %>%
  mutate(period = factor(
    period,
    levels = c(
      "extreme_sum1",
      "extreme_sum5",
      "extreme_sum10",
      "extreme_sum15",
      "extreme_sum20"
    ),
    labels =
      c(
        TeX("$prop_{i,12}^{(12)}$"),
        TeX("$prop_{i,60}^{(12)}$"),
        TeX("$prop_{i,120}^{(12)}$"),
        TeX("$prop_{i,180}^{(12)}$"),
        TeX("$prop_{i,240}^{(12)}$")
      )
  )) %>%
  filter(!is.na(extreme_sum)) %>%
  ggplot(aes(x = Datum, y = extreme_sum)) +
  geom_line() +
  labs(x = "Year", y = "Number of spei12 less then -2 \n") +
  facet_wrap(
    ~ period,
    ncol = 1,
    strip.position = "left",
    labeller = label_parsed
  ) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    panel.spacing = unit(0.1, "lines"),
    axis.title.y = element_blank()
  )

# remove everything from the global environment except the plots and wasserstein_distances
rm(list = setdiff(
  ls(),
  c("climate_data",
    "modelfitting_plot",
    "wasserstein_distances",
    "spei_plot",
    "extreme_droughts_plot",
    "extreme_droughts_time_scales_plot",
    "extreme_droughts_periods_plot"
  )
))