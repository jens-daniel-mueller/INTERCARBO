#### load required packages ####

library(data.table)
library(ggplot2)
library(tidyverse)
library(lubridate)



#### Define function to calculate pHT from spectrophotometric measurements ####
#### according to Mueller and Rehder (2018)

pHT.Mueller <- function(Sal, Tem, Rspec){
  
  #first set of coefficients defines pK2e2 = f(Sal, Tem)
  
  1.08071477e+03                      -
    1.35394946e-01  *Sal^0.5            -   
    1.98063716e+02  *Sal^1.5            +
    6.31924397e+01  *Sal^2              -
    5.18141866e+00  *Sal^2.5            -
    2.66457425e+04  *Tem^-1             +
    5.08796578e+03  *Sal^1.5 * Tem^-1   -
    1.62454827e+03  *Sal^2 * Tem^-1     +
    1.33276788e+02  *Sal^2.5 * Tem^-1   -
    1.89671212e+02  *log(Tem)           +
    3.49038762e+01  *Sal^1.5 * log(Tem) -
    1.11336508e+01  *Sal^2 * log(Tem)   +
    9.12761930e-01  *Sal^2.5 * log(Tem) +
    3.27430677e-01  *Tem              -
    7.51448528e-04  *Sal^0.5 * Tem      +
    3.94838229e-04  *Sal * Tem          -
    6.00237876e-02  *Sal^1.5 * Tem      +
    1.90997693e-02  *Sal^2 * Tem        -
    1.56396488e-03  *Sal^2.5 * Tem      +
    
    #second set of coefficients includes the definition of mCP absorptivity ratios e1 and e3/e3
    #as determined by Liu et al. (2011) and defines the log-term calculation 
    
    log10(
      (Rspec -
         (-0.007762 + 4.5174e-5*Tem)) /
        (1 - (Rspec *  (- 0.020813 + 2.60262e-4*Tem + 1.0436e-4*(Sal-35))))
    )
}

#### Define function to calculate target pHT of TRIS buffer solutions ####

pHT.TRIS <- function(Sal, Tem, b_Tris){
  -327.3307 -
    2.400270 * Sal +
    8.124630e-2 * Sal^2 -
    9.635344e-4 * Sal^3 -
    
    9.103207e-2 * Tem -
    1.963311e-3 * Sal * Tem +
    6.430229e-5 * Sal^2 * Tem -
    7.510992e-7 * Sal^3 * Tem +
    
    56.92797 * log(Tem) +
    5.235889e-1 * Sal * log(Tem) -
    1.7602e-2 * Sal^2 * log(Tem) +
    2.082387e-4 * Sal^3 * log(Tem) +
    
    11382.97 * (1/Tem) -
    
    2.417045 * b_Tris +
    7.645221e-2 * b_Tris * Sal +
    1.122392e-2 * b_Tris * Tem -
    3.248381e-4 * b_Tris * Sal * Tem -
    
    4.161537 * b_Tris^2 +
    6.143395e-2 * b_Tris^2 * Sal
}


# Test value: pHT = 8.0703 at Sal = 20, Tem = 298.15, and B_Tris = 0.04
pHT.TRIS(20, 298.15, 0.04)


#### Load raw data files, select and rename columns #### 

setwd("C:/Mueller_Jens_Data/Projects/181122_TNA_Oslo_intercomparison_pH_pCO2/data/Carter_spec_pH")
df <- read.delim("181125_Oslo_Mueller_Datafile",  header=FALSE)

df <- data.table (df[,c(seq(3,7,1), seq(9,13,1), 15)])
names(df) <- 
  c("date", "time", "dye", "sample", "Rep", "A5", "A4", "A7", "Ai", "Tem", "Sal")

#### Correct date of samples that were measured the next day after sampling by -1 ####

df[sample=="T10S35PC800-2038"]$date <- "11/23/2018"
df[sample=="T10S20PC800-2120"]$date <- "11/24/2018"
df[sample=="T10S05PC200-2145"]$date <- "11/25/2018"


#### Remove measurements done on test samples ####

df <- df[sample != "T10S35PC400-1509"]


#### Subset data frame and calculate R value ####

df <- df[Rep != 0]
df$Rspec <- (df$A5-df$A7) / (df$A4-df$A7)

#### Calculate sample pHT values ####

df$pHT.Mueller <- pHT.Mueller(df$Sal, df$Tem, df$Rspec)

#### Subset data set for buffer solutions (TRIS) and subsample from test tank (SW) ####

df$solution <- substr(df$sample, 1, 3)
TRIS <- df[solution == "TRI"]
SW <- df[solution != "TRI"]

rm(df)

#### Evaluate measurements performed on PTB TRIS buffer solutions at S = 5, 20, and 35 ####

#### calculate target pHT values

TRIS$pHT.target <- pHT.TRIS(TRIS$Sal, TRIS$Tem, 0.04)

#### Subset TRIS data frame for PTB and IOW buffer solutions

TRIS$Source <- substr(TRIS$sample, 6, 8)
TRIS <- TRIS[Source != "PTB" | Ai < 0.7] #remove one measurement with too high dye concentration

TRIS.PTB <- TRIS[Source == "PTB"]
TRIS.IOW <- TRIS[Source == "IOW"]

rm(TRIS)

#### Plots

TRIS.PTB %>% 
  filter(Tem>297) %>% 
  ggplot()+
  geom_point(aes(Ai, pHT.target, col=date, shape="target"))+
  geom_point(aes(Ai, pHT.Mueller, col=date, shape="observed"))+
  facet_wrap(~Sal, scales = "free_y")+
  xlim(0,0.6)+
  scale_y_continuous(breaks = seq(5,10,0.005))

TRIS.PTB %>% 
ggplot()+
  geom_point(aes(Ai, pHT.Mueller - pHT.target, col=date))+
  facet_wrap(~Sal, scales = "free_y")+
  xlim(0,0.6)+
  scale_y_continuous(breaks = seq(-5,10,0.002))

TRIS.IOW %>% 
ggplot()+
  geom_point(aes(Ai, pHT.target, col=date, shape="target"))+
  geom_point(aes(Ai, pHT.Mueller, col=date, shape="obs"))+
  facet_grid(round(Tem, 0)~Sal, scales = "free_y")+
  xlim(0,0.6)+
  scale_y_continuous(breaks = seq(5,10,0.005))


rm(TRIS.IOW, TRIS.PTB)


#### Analyse reference pH measurements of samples from test tanks ####

#### Create individual columns for target Temperature, Salinity, and pCO2

SW <-
data.table(
  SW %>% 
  mutate(Tem.target = substr(sample, 2, 3), 
         Sal.target = substr(sample, 5, 6), 
         pCO2.target = substr(sample, 9, 11),
         time.sample = substr(sample, 13, 16)))

SW[time.sample == "NaOH"]$time.sample <- substr( SW[time.sample == "NaOH"]$sample, 18, 21)
SW$solution <- NULL  

#### Clean data set from errornous measurements and restrict to measurements at 25oC

SW.25 <-  
SW %>% 
  filter(Tem == 298.15,                       #remove comparison measurement at in-situ temperature
         sample != "T20S35pC200-1310",        #remove measurements of T20S35pC200 before NaOH addition
         sample != "T10S35pC400-1509",        #omit inital test measurements of T10S35pC400-1509
         Rep < 4,                             #remove last replicate potentially affected by gas exchange 
         sample != "T20S20pC200-1050" | Rep != 1 | dye != "p")



#### transform time and data columns to POSIXct

SW$time.sample <- mdy_hm(paste(SW$date, SW$time.sample ))
SW.25$time.sample <- mdy_hm(paste(SW.25$date, SW.25$time.sample ))


#### Extrapolate measurements to zero dye concentration
  

SW.25.ex <-
SW.25 %>% 
  group_by(Tem.target, Sal.target, pCO2.target, time.sample, dye) %>% 
  summarise(pHT.ex = lm(pHT.Mueller ~ Ai)$coefficients[[1]],
            slope  = lm(pHT.Mueller ~ Ai)$coefficients[[2]],
            SD.residuals = sd(lm(pHT.Mueller ~ Ai)$residuals) )  %>% 
  ungroup()



SW.25 %>%   
  ggplot()+
  geom_smooth(method = "lm", aes(Ai, pHT.Mueller, linetype=dye), 
              se=FALSE, col = "grey50", fullrange = TRUE)+
  geom_point(aes(Ai, pHT.Mueller, col=dye, shape=as.factor(Rep) ))+
  geom_point(data=SW.25.ex, aes(0, pHT.ex, col=dye, shape="ex" ))+
  facet_wrap(~Tem.target*Sal.target*pCO2.target*time.sample, scales = "free_y", ncol =3)+
  scale_color_brewer(palette = "Set1")+
  scale_y_continuous(breaks = seq(5,10,0.005)) 


SW.25.ex %>% 
ggplot(aes(pHT.ex, slope, col=dye))+
  geom_hline(yintercept = 0)+
  geom_point()+
  geom_errorbar(aes (ymin = slope - SD.residuals, ymax= slope + SD.residuals))+
  facet_wrap(~Sal.target)



SW.25.ex.save <-
  SW.25.ex %>% 
  filter(dye == "p") %>% 
  select(-c(dye, slope, SD.residuals)) %>% 
  arrange(time.sample)


#### Safe final data set ####


setwd("C:/Mueller_Jens_Data/Projects/181122_TNA_Oslo_intercomparison_pH_pCO2/data/_Finalized_datasets")
write.csv(SW.25.ex.save, "IOW_Carter_pH.csv", row.names = FALSE)
  
rm(list=ls())

