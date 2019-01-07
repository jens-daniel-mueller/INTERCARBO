library(tidyverse)
library(lubridate)
library(data.table)
library(plotly)


####

setwd("C:/Mueller_Jens_Data/Projects/181122_TNA_Oslo_intercomparison_pH_pCO2/data/HydroFIA_pH/IOW_data")

HF1 <- read.delim("1_PH-0218-001-data.txt", sep=",", skip = 2)[,c(1,3,6,7,8,9)]
HF2 <- read.delim("2_PH-0515-001-data.txt", sep=",", skip = 2)[,c(1,3,6,7,8,9)]
HF3 <- read.delim("3_PH-1017-001-data.txt", sep=",", skip = 2)[,c(1,3,6,7,8,9)]

HF1$instrument <- "1_PH-0218-001"
HF2$instrument <- "2_PH-0515-001"
HF3$instrument <- "3_PH-1017-001"

HF <- data.table( rbind(HF1, HF2, HF3) )
rm(HF1, HF2, HF3)

names(HF) <- c("date", "sample", "Sal", "pH.Mosley.25", "pH.fitpoints", "pH.error", "instrument")

HF$date <- ymd_hms(HF$date)

HF <-
HF %>% 
mutate(Tem.target = substr(sample, 2, 3), 
       Sal.target = substr(sample, 5, 6), 
       pCO2.target = substr(sample, 9, 11))


ggplot()+
  geom_path(data = HF, aes(date, pH.Mosley.25, col=instrument))+
  geom_point(data = SW, aes(time.sample, pHT.Mueller))+
  scale_x_datetime(date_labels = "%d.%b, %H:%M")+
  facet_wrap(~Tem.target*Sal.target*pCO2.target, scales = "free", ncol =3)
