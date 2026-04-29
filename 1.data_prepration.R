#rm(list = ls(all=T))

library(haven)
library(tidyverse)
library(descr)

# loading the data -----
uk<- readRDS("data/mhc_can_data.rds")

## Re-coding variables ------
uk$menopause<- factor(uk$menopause,
                           labels = c("no", "yes", "Not applicable"))

uk$blm.t2d<- factor(uk$blm.t2d, 
                         labels = c("no", "yes"))

uk<- uk %>% mutate(smok= case_when(smoking=="Current"~ "Current", 
                                             smoking=="Never"~ "Never",
                                             smoking=="Previous" ~"Previous" 
))

uk<- uk %>% mutate(ethn= case_when(ethnicity=="White"~"White",
                                             ethnicity=="Black"~"Black",
                                             ethnicity=="Chinese"~"Chinese",
                                             ethnicity=="Mixed"~"Mixed",
                                             ethnicity=="South Asian"~"South Asian", 
                                             ethnicity=="any other" ~"anyother"))

uk$ethn<- factor(uk$ethn, levels = c("White", "Black","Chinese", "Mixed",
                                               "South Asian","anyother" ))


uk<-uk %>% mutate(sex= case_when(sex== "M" ~"M", 
                                           sex== "F" ~ "F", 
                                           NA~NA))

uk$sex<- as.factor(uk$sex)


uk$ffq.pmi <- factor(uk$ffq.pmi, labels = c("Never", "Less than once a week", 
                                            "Once a week", "2-5 times a week", 
                                            "5-6 times a week", "once or more daily" ))

#Creating a new MHCs varible combining all conditions
uk <-uk %>% mutate(mhc.e= 
                     if_else(dep.e==1 | anx.e ==1 | bpl.e==1 |
                               srd.pts.e == 1 | sci.e== 1, 1,0))

#Filtering to remove incident cancer before the baseline 
uk <- uk %>% filter(can.t>=0)


#rm(list = setdiff(ls(), "uk"))
