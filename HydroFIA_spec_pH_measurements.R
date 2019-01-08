#### load required packages ####
library(tidyverse)
library(lubridate)
library(data.table)
library(seacarb)

#### load processed Carter pH data file ####

setwd("C:/Mueller_Jens_Data/Projects/181122_TNA_Oslo_intercomparison_pH_pCO2/data/_Finalized_datasets")
SW.25.ex.save <- read_csv("IOW_Carter_pH.csv")

#### load HydroFIA data files ####

setwd("C:/Mueller_Jens_Data/Projects/181122_TNA_Oslo_intercomparison_pH_pCO2/data/HydroFIA_pH/IOW_data")

# # original HydroFIA data files with pH calculated according to outdated Mosley et al (2004)
# HF1 <- read.delim("1_PH-0218-001-data.txt", sep=",", skip = 2)[,c(1,3,6,7,8,9)]
# HF2 <- read.delim("2_PH-0515-001-data.txt", sep=",", skip = 2)[,c(1,3,6,7,8,9)]
# HF3 <- read.delim("3_PH-1017-001-data.txt", sep=",", skip = 2)[,c(1,3,6,7,8,9)]

# re-processed HydroFIA data files with pH calculated according to updated Mueller and Rehder (2018)
HF1 <- read_csv("1_PH-0218-001-pH-recalculated.csv") %>% 
  mutate(instrument = "1_PH-0218-001")
HF2 <- read_csv("2_PH-0515-001-pH-recalculated.csv") %>% 
  mutate(instrument = "2_PH-0515-001")
HF3 <- read_csv("3_PH-1017-001-pH-recalculated.csv") %>% 
  mutate(instrument = "3_PH-1017-001")

HF <- bind_rows(HF1, HF2, HF3)
rm(HF1, HF2, HF3)

# names(HF) <- c("date", "sample", "Sal", "Tem", "counter",
#                "pH.Mueller", "pH.Mueller.error", "pH.Mosley", "pH.Mosley.error", "instrument")


#### Create individual columns for target temperature, salinity, and pCO2 ####

HF <-
HF %>% 
mutate(Tem.target = substr(sampleName, 2, 3), 
       Sal.target = substr(sampleName, 5, 6), 
       pCO2.target = substr(sampleName, 9, 11))

#### Subset measurements performed on tanks ####

HF.sw <-
  HF %>%
  filter(Tem.target != "RI")

#### Assign Alkalinity values and correct pH readings to temperature = 25oC ####

HF.sw <-
  HF.sw %>% 
  mutate(AT = (salinity/35)*2300e-6,
         pH.Mueller.25 = pHinsi(pH = ph_mueller, ALK = AT, Tinsi = 25, Tlab = temperature))


#### Remove initial measurements before stabile results are achieved ####

HF.sw.stabil <-
  HF.sw %>%
  filter(timeStamp > dmy_hm("25/11/2018 12:00") | sampleName != "T10S05PC400") %>% 
  filter(timeStamp > dmy_hm("24/11/2018 14:30") | sampleName != "T10S20PC400") %>% 
  filter(timeStamp > dmy_hm("24/11/2018 10:30") | sampleName != "T10S35PC200") %>% 
  filter(timeStamp > dmy_hm("26/11/2018 16:15") | sampleName != "T20S05PC200")


#### Plot results ####

setwd("C:/Mueller_Jens_Data/Projects/181122_TNA_Oslo_intercomparison_pH_pCO2/plots")

 ggplot()+
  geom_path(data = HF.sw.stabil, aes(timeStamp, pH.Mueller.25, col=instrument))+
  geom_point(data = SW.25.ex.save, aes(time.sample, pHT.ex, col="Carter"))+
  scale_x_datetime(date_labels = "%d.%b, %H:%M")+
  scale_color_brewer(palette = "Set1", name="Instrument")+
  facet_wrap(~Tem.target*Sal.target*pCO2.target,
             scales = "free", ncol =3, labeller = label_both)+
   labs(x="Date", y="pHT @ 25oC")+
   theme_bw()

ggsave("IOW_Carter_and_HydroFIA_pH_facetted.jpg", height = 300, width = 500, units = "mm")


ggplot()+
  geom_path(data = HF.sw, aes(timeStamp, pH.Mueller.25, group=instrument), col="grey80")+
  geom_path(data = HF.sw.stabil, aes(timeStamp, pH.Mueller.25,
                                     col=instrument, group=interaction(sampleName, instrument)))+
  geom_point(data = SW.25.ex.save, aes(time.sample, pHT.ex, col="Carter"))+
  scale_x_datetime(date_labels = "%d.%b, %H:%M")+
  scale_color_brewer(palette = "Set1", name="Instrument")+
  labs(x="Date", y="pHT @ 25oC")+
  theme_bw()

ggsave("IOW_Carter_and_HydroFIA_pH_timeseries.jpg", height = 210, width = 297, units = "mm")


# HF.sw %>% 
#   ggplot()+
#   geom_histogram(aes(temperature))
