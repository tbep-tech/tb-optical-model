# enable repos
options(repos = c(
  tbeptech = 'https://tbep-tech.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))

# install tbeptools
install.packages('tbeptools')

# load tbeptools
library(tbeptools)
library(tidyverse)

setwd("./data-example/")
xlsx <- "RWMSpreadsheet_ThroughCurrentReportMonth.xlsx"
epcdata <- read_importwq(xlsx, download_latest = TRUE)

test_batch <- epcdata %>%
              dplyr::filter(yr==2026) %>%
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
plot_seasonal_KdPAR(unlist(results))
plot_Kd_spectrum(unlist(results))         #Screwy plotresults here?
plot_absorption(Unlist(results))          #Plots not working here?
