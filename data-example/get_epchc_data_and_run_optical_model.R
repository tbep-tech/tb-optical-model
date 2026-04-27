# enable repos
options(repos = c(
  tbeptech = 'https://tbep-tech.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))

# install tbeptools
install.packages(c('tbeptools', 'tidyverse', 'here'))



# Load R libraries
library(tbeptools)
library(tidyverse)
library(here)

source('./R/original_recode.R')

setwd("./data-example/")
xlsx <- "RWMSpreadsheet_ThroughCurrentReportMonth.xlsx"

epcall <- read_importepc(xlsx, download_latest = TRUE)

epcdata <- read_importwq(xlsx, download_latest = TRUE)

test_batch <- epcdata %>%
              dplyr::filter(yr==2025) %>%
              mutate(Station_ID = epchc_station,
                     Year = yr,
                     Month = mo,
                     Day = day(SampleTime),
                     DecTime = hour(SampleTime)+(minute(SampleTime)/60),
                     z = Total_Depth_m,
                     CDOM440 = as.numeric(Color_345_F45_PCU),
                     CHLA = as.numeric(chla),
                     NTU = as.numeric(`Turbidity_JTU-NTU`)) %>%
              select(Station_ID, Year, Month, Day, DecTime, z, CDOM440, CHLA, NTU)

results <- batch_run(test_batch)
plot_seasonal_KdPAR(results)
plot_Kd_spectrum(results[215, ])
plot_absorption(results[215,7], results[215,8], results[215,9])
sensitivity_analysis(results[215,7], results[215,8], results[215,9], results[215,10])
