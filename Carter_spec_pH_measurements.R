library(data.table)
library(ggplot2)
library(tidyverse)
library(lubridate)


#### load raw data files, select and rename columns #### 

setwd("C:/Mueller_Jens_Data/Projects/181122_TNA_Oslo_intercomparison_pH_pCO2/data/Carter_spec_pH")
df <- read.delim("181125_Oslo_Mueller_Datafile",  header=FALSE)

df <- data.table (df[,c(seq(3,7,1), seq(9,13,1), 15)])
names(df) <- 
  c("date", "time", "dye", "sample", "Rep", "A5", "A4", "A7", "Ai", "Tem", "Sal")


#### Subset data frame and calculate R value ####

df <- df[Rep != 0]
df$Rspec <- (df$A5-df$A7) / (df$A4-df$A7)


#### Definitions of equations to calculate pH from measured ####
#### Salinity (Sal), Temperature in K (Tem) and mCP absorption ratio (Rspec = A434/A578)

#### Mosley et al. 2004

pHT.Mosley <- function(Sal, Tem, Rspec,
                       e1 = 0.00691, e2 = 2.222, e3 = 0.1331) {
  (1245.69/Tem) + 4.4572353 - (0.3238*(Sal^0.5)) + (0.0807*Sal) - (0.01157*(Sal^1.5)) + (0.000694*(Sal^2)) +
    log10( (Rspec-e1) / (e2-Rspec*e3) ) 
}


#### Mueller and Rehder (2018)

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


#### Calculate sample pHT value according to two mCP parameterizations ####

df$pHT.Mueller <- pHT.Mueller(df$Sal, df$Tem, df$Rspec)
df$pHT.Mosley <- pHT.Mosley(df$Sal, df$Tem, df$Rspec)


#### subset data set for buffer solutions (TRIS) and subsample from test tank (SW) ####

df$solution <- substr(df$sample, 1, 3)

TRIS <- df[solution == "TRI"]
SW <- df[solution != "TRI"]

rm(df)



#### Evaluate measurements performed on PTB TRIS buffer solutions at S = 5, 20, and 35 ####

#### pHT model of TRIS buffer solutions to calculate target values ####

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


#### Subset TRIS data frame for high quality solutions provided by Frank Bastkowski (PTB) ####

TRIS$Source <- substr(TRIS$sample, 6, 8)
TRIS <- TRIS[Source != "PTB" | Ai < 0.7] #remove one measurement with too high dye concentration

TRIS$pHT.true <- pHT.TRIS(TRIS$Sal, TRIS$Tem, 0.04)

TRIS[Source == "PTB"] %>% 
  ggplot()+
  geom_point(aes(Ai, pHT.true, col=date, shape="true"))+
  geom_point(aes(Ai, pHT.Mueller, col=date, shape="obs"))+
  facet_wrap(~Sal, scales = "free_y")+
  xlim(0,0.6)+
  scale_y_continuous(breaks = seq(5,10,0.005))

TRIS[Source == "PTB"] %>% 
  ggplot()+
  geom_point(aes(Ai, pHT.Mueller - pHT.true, col=date))+
  facet_wrap(~Sal, scales = "free_y")+
  xlim(0,0.6)+
  scale_y_continuous(breaks = seq(-5,10,0.002))


#### Subset TRIS data frame for IOW homemade solutions ####

TRIS[Source == "IOW"] %>% 
  ggplot()+
  geom_point(aes(Ai, pHT.true, col=date, shape="true"))+
  geom_point(aes(Ai, pHT.Mueller, col=date, shape="obs"))+
  facet_grid(round(Tem, 0)~Sal, scales = "free_y")+
  xlim(0,0.6)+
  scale_y_continuous(breaks = seq(5,10,0.005))



SW <-
data.table(
  SW %>% 
  mutate(Tem.target = substr(sample, 2, 3), 
         Sal.target = substr(sample, 5, 6), 
         pCO2.target = substr(sample, 9, 11),
         time.sample = substr(sample, 13, 16)))

SW[time.sample == "NaOH"]$time.sample <- substr( SW[time.sample == "NaOH"]$sample, 18, 21)
SW$solution <- NULL  

SW.25 <-  
SW %>% 
  filter(Tem == 298.15,                       #remove comparison measurement at in-situ temperature
         sample != "T20S35pC200-1310",        #remove measurements of T20S35pC200 before NaOH addition
         sample != "T10S35pC400-1509",        #omit inital test measurements of T10S35pC400-1509
         Rep < 4,                             #remove last replicate potentially affected by gas exchange 
         sample != "T20S20pC200-1050" | Rep != 1 | dye != "p")
  
SW.25 %>%   
  ggplot()+
  geom_smooth(method = "lm", aes(Ai, pHT.Mueller, linetype=dye), 
              se=FALSE, col = "grey50", fullrange = TRUE)+
  geom_point(aes(Ai, pHT.Mueller, shape=dye, col=as.factor(Rep) ))+
  facet_wrap(~Tem.target*Sal.target*pCO2.target, scales = "free_y", ncol =3)+
  xlim(0,0.8)+
  scale_color_brewer(palette = "Set1")+
  scale_y_continuous(breaks = seq(5,10,0.005)) 

extrapol <-
SW.25 %>% 
  group_by(Tem.target, Sal.target, pCO2.target, dye) %>% 
  summarise(pHT.ex := lm(pHT.Mueller ~ Ai)$coefficients[[1]],
            slope  := lm(pHT.Mueller ~ Ai)$coefficients[[2]],
            SD.residuals := sd(lm(pHT.Mueller ~ Ai)$residuals) )  %>% 
  ungroup()


ggplot(extrapol, aes(pHT.ex, slope, col=dye))+
  geom_hline(yintercept = 0)+
  geom_point()+
  geom_errorbar(aes (ymin = slope - SD.residuals, ymax= slope + SD.residuals))+
  facet_wrap(~Sal.target)



SW.25 %>%   
  ggplot()+
  geom_smooth(method = "lm", aes(Ai, pHT.Mueller, linetype=dye), 
              se=FALSE, col = "grey50", fullrange = TRUE)+
  geom_point(aes(Ai, pHT.Mueller, shape=dye, col=as.factor(Rep) ))+
  facet_wrap(~Tem.target*Sal.target*pCO2.target, scales = "free_y", ncol =3)+
  scale_color_brewer(palette = "Set1")+
  scale_y_continuous(breaks = seq(5,10,0.005))+
  geom_point(data=extrapol, aes(0, pHT.ex))


SW$time.sample <- mdy_hm(paste(SW$date, SW$time.sample ))















