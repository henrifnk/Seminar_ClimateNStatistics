library(tidyverse)
library(SPEI)
library(zoo)

# Load the temperature data
for (i in 1:12) {
  file_name <- paste0("work/02-spei/data/regional_averages_tm_", sprintf("%02d", i), ".txt")
  data <- read_delim(file_name, delim = ";", escape_double = FALSE, trim_ws = TRUE, skip = 1)
  data <- data %>% dplyr::select(Jahr, Monat, Bayern)
  assign(paste0("data_tm_", sprintf("%02d", i)), data)
}

for (i in 1:12) {
  assign(
    paste0("data_tm_", sprintf("%02d", i)),
    transmute(
      get(
        paste0("data_tm_", sprintf("%02d", i))),
      temperature = as.numeric(Bayern),
      Jahr = as.numeric(Jahr),
      Monat = as.numeric(Monat)))
}

data_tm <- bind_rows(data_tm_01, data_tm_02, data_tm_03, data_tm_04,
                         data_tm_05, data_tm_06, data_tm_07, data_tm_08,
                         data_tm_09, data_tm_10, data_tm_11, data_tm_12) %>%
  arrange(Jahr, Monat)

rm(list = ls(pattern = "data_tm_"))

# Load the precipitation data
for (i in 1:12) {
  file_name <- paste0("work/02-spei/data/regional_averages_rr_", sprintf("%02d", i), ".txt")
  data <- read_delim(file_name, delim = ";", escape_double = FALSE, trim_ws = TRUE, skip = 1)
  data <- data %>% dplyr::select(Jahr, Monat, Bayern)
  assign(paste0("data_rr_", sprintf("%02d", i)), data)
}

for (i in 1:12) {
  assign(
    paste0("data_rr_", sprintf("%02d", i)),
    transmute(
      get(
        paste0("data_rr_", sprintf("%02d", i))),
      precipitation = as.numeric(Bayern),
      Jahr = as.numeric(Jahr),
      Monat = as.numeric(Monat)))
}

data_rr <- bind_rows(data_rr_01, data_rr_02, data_rr_03, data_rr_04,
                         data_rr_05, data_rr_06, data_rr_07, data_rr_08,
                         data_rr_09, data_rr_10, data_rr_11, data_rr_12) %>%
  arrange(Jahr, Monat)

rm(list = ls(pattern = "data_rr_"))

# Merge the temperature and precipitation data
climate_data <- left_join(data_tm, data_rr, by = c("Jahr", "Monat"))

# Latitude of the region
lat_bayern <- 48.1351

# Create a date column
climate_data$Datum <- as.Date(paste0(climate_data$Jahr, "-", sprintf("%02d", climate_data$Monat), "-01"))

# Create a time series object for the temperature
climate_data$temperature <- ts(climate_data$temperature, frequency = 12, start = c(1881, 1))

# Calculate PET using the thornthwaite method implemented in the SPEI package
climate_data$PET <- thornthwaite(climate_data$temperature, lat = lat_bayern)

# D is the difference between precipitation and PET (D_i = P_i - PET_i)
climate_data$D <- climate_data$precipitation - climate_data$PET

# define the reference period
ref_start <-c(1980, 1)
ref_end <- c(2010, 12)

# define the time scales
time_scales <- c(1, 3, 6, 12, 18, 24)

# calculate the SPEI values at different time scales
spei_data <- map(time_scales, ~ spei(climate_data$D, scale = .x, ref.start = ref_start, ref.end = ref_end)$fitted)
climate_data <- climate_data <- climate_data %>%
  mutate(!!!setNames(spei_data, paste0("spei", time_scales)))

# The D_i values are aggregated at different time scales
climate_data <- climate_data %>%
  mutate(!!!setNames(map(time_scales, ~ rollapply(climate_data$D, .x, sum, fill = NA, align = "right")), paste0("D_accum_", time_scales)))

climate_data <- climate_data %>%
  select(Datum, Jahr, Monat, everything())

# export the data
saveRDS(climate_data, "work/02-spei/data/climate_data.rds")