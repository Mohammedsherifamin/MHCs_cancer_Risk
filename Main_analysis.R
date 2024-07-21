## Main analysis
# Load packages and data

library(haven)
library(tidyverse)
library(descr)
library(forcats)
library(foreign)
library(gtsummary)
library(lattice)
library(mosaic)
library(broom)
library(Hmisc)
library(lmtest)
library(survival)
library(flextable)

ds<- read_csv("ukb_project_2023.csv")
extra<- read_csv("extra.csv")

## Data preparation----
#Merging the two datasets by participant identification code 
ds <- merge(ds, extra, by = "eid", all.x = TRUE)

#Creating a new MHCs variable combining individual conditions
ds<- ds %>% mutate(mhc.e= 
                       if_else(dep.e==1 | anx.e ==1 | bpl.e==1 |
                                 srd.pts.e == 1 | sci.e== 1, 1,0), 
                     mhc.t = 
                       pmin(dep.t, anx.t, bpl.t, srd.pts.t, sci.t, na.rm = TRUE))

## Re-coding variables to prepare for regression
ds$menopause<- factor(ds$menopause, 
                       labels = c("no", "yes", "Not applicable"))

ds$blm.t2d<- factor(ds$blm.t2d, 
                     labels = c("no", "yes"))

ds<- ds %>% mutate(smok= case_when(smoking=="Current"~ "Current", 
                                     smoking=="Never"~ "Never",
                                     smoking=="Previous" ~"Previous" 
))

ds<- ds %>% mutate(ethn= case_when(ethnicity=="White"~"White",
                                     ethnicity=="Black"~"Black",
                                     ethnicity=="Chinese"~"Chinese",
                                     ethnicity=="Mixed"~"Mixed",
                                     ethnicity=="South Asian"~"South Asian", 
                                     ethnicity=="any other" ~"anyother"))

ds$ethn<- factor(ds$ethn, levels = c("White", "Black","Chinese", "Mixed",
                                     "South Asian","anyother" ))


ds<-ds %>% mutate(sex= case_when(sex== "M" ~"M", 
                                   sex== "F" ~ "F", 
                                   NA~NA))

ds$sex<- as.factor(ds$sex)

uk<- uk %>% mutate(
  mhc.e_new = case_when(srd.pts.e==1~ "pts",
                        dep.e==1~ "dep", 
                        anx.e==1~ "anx",
                        bpl.e==1~ "bpl", 
                        sci.e==1~ "sci"))

uk$mhc.e_new <- factor(uk$mhc.e_new, levels = c("dep", "anx","bpl", "sci",
                                                "pts"))

## Check the amount of missing data 
missing_counts <- colSums(is.na(ds))   

## Due to the small number of missing data, I decided to drop missing data and preform a 
#complete cases analysis
ds <- ds %>%  drop_na()

#Filtering cancer to above baseline 
ds <- ds %>% filter(can.t>0)

## Table summary of  baseline characteristics -----
base_char<-  c("age", "sex", "ethn", "deprivation" , 
               "smok", "alc.ut" , "slp.durat",   "met.tot",  "ffq.pmi",
               "bmi" , "blm.t2d",   "sbp"    , "menopause", "can.t")

tbl.sum <- tbl_summary(ds, missing = "no",
                       include = base_char , 
                       statistic = list( all_continuous() ~ "{mean}, ({sd})", 
                                         c("can.t", "alc.ut") ~ "{median} ({p25}, {p75})"),
                       type = list(c(ffq.pmi, slp.durat) ~ "continuous")) 

tbl.sum %>% as_flex_table() %>% save_as_docx(path = "tab/summary_overall.docx")

#Stratifying baseline characteristics by sex
tbl.sum.sex <- tbl_summary(ds,  by = sex, 
                           missing = "no", 
                           include = base_char , 
                           statistic = list( all_continuous() ~ "{mean}, ({sd})", 
                                             c("can.t", "alc.ut") ~ "{median} ({p25}, {p75})"),
                           type = list(c(ffq.pmi, slp.durat) ~ "continuous")) %>% 
  add_overall() 

# Testing for statistical differences between groups
tbl.sum.sex <- tbl.sum.sex %>% 
  add_p(test = list(all_continuous() ~ "t.test", c(can.t, alc.ut) ~ "wilcox.test"))

tbl.sum.sex %>% as_flex_table() %>% save_as_docx(path = "tab/summary_sex.docx")


#Stratifying baseline characteristics by by MHCs diagnosis
tbl.sum.mhc <- tbl_summary(ds,  by = mhc.e, 
                           missing = "no", 
                           include = base_char , 
                           statistic = list( 
                             all_continuous() ~ "{mean}, ({sd})", 
                             c("can.t", "alc.ut") ~ "{median} ({p25}, {p75})"),
                           type = list(c(ffq.pmi, slp.durat) ~ "continuous")) %>% 
  add_overall() 

tbl.sum.mhc <- tbl.sum.mhc %>% 
  add_p(test = list(all_continuous() ~ "t.test", 
                    c(can.t , alc.ut) ~ "wilcox.test"))

print(tbl.sum.mhc)

tbl.sum.mhc %>% as_flex_table() %>% 
  save_as_docx(path = "tab/summary_mhc.docx")

#Table summary for PHQ-4 scores 
# table for PHQ-4
tbl.phq <- tbl_summary(ds,  by = mhc.e_new, 
                       missing = "no", 
                       include = c("phq") , 
                       statistic = list(c("phq") ~ "{median} ({p25}, {p75})")) %>% 
  add_overall() 


tbl.phq <- tbl.phq %>% 
  add_p()


tbl.phq %>% as_flex_table() %>% save_as_docx(path = "tab/phq.docx")

## Cox Regression Models with overall cancer ------
#Depression and overall cancer risk-----

#Creat a new data set
ds <- ds %>% 
  # a new variable to ensure a 1 year land mark analysis and that post baseline depression cases
  # are only those with depression diagnosis before cancer diagnosis 
  mutate(dep.e_new = if_else(dep.t - can.t < -1 , 1, 0), 
         can.t_new = if_else(dep.e_new==1 & dep.t > 0, can.t - dep.t, can.t))

ds$dep.e_new<- as.factor(ds$dep.e_new)

#Adjusted for Demographics "model 1"
depcox1<- coxph(Surv(can.t_new, can.e)~ dep.e_new +
                  age +sex + ethn + deprivation, 
                data = ds)

#Adjusted for Socidemographics + lifestyle "model 2"
depcox2<- coxph(Surv(can.t_new, can.e)~ dep.e_new +
                  age +sex + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi ,
                data = ds)

#FullyAdjusted "model 3"
depcox3<- coxph(Surv(can.t_new, can.e) ~ 
                  dep.e_new + age +sex + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp , data=ds)

#Stratified by sex

fdepcox<- coxph(Surv(can.t_new, can.e) ~ dep.e_new +
                  age +sex + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp + menopause , 
                data= subset(ds, sex=="F")) 


mdepcox<- coxph(Surv(can.t_new, can.e) ~ dep.e_new +
                  age +sex + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp , data= subset(ds, sex=="M")) 

## Anxiety and overall cancer risk ----
ds <- ds %>% 
  mutate(anx.e_new = if_else(anx.t - can.t < -1, 1, 0),
         can.t_new = if_else(anx.e_new==1 & anx.t > 0, can.t - anx.t, can.t))

ds$anx.e_new<- as.factor(ds$anx.e_new)

#model 1
anxcox1<- coxph(Surv(can.t_new, can.e) ~ anx.e_new+ 
                  age +sex + ethn + deprivation, 
                data=ds) 

#model 2 
anxcox2<- coxph(Surv(can.t_new, can.e) ~ anx.e_new + 
                  age +sex + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                data=ds) 

#model 3
anxcox3<- coxph(Surv(can.t_new, can.e) ~ anx.e_new + 
                  age +sex + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp  , data=ds) 

## Bipolar Disorder and overall cancer risk -----
ds<- ds %>%  
  filter(can.t>0) %>% 
  mutate(bpl.e_new = if_else(bpl.t - can.t < -1, 1, 0),
         can.t_new = if_else(bpl.e_new==1 & bpl.t >0, can.t - bpl.t, can.t))

ds$bpl.e_new<- factor(ds$bpl.e_new)

#model 1
bplcox1<- coxph(Surv(can.t_new, can.e)~  bpl.e_new  + age +sex + ethn + 
                  deprivation , data=ds)

#model 2
bplcox2<- coxph(Surv(can.t_new, can.e)~  bpl.e_new  +
                  age +sex + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                data=ds)

#model 3
bplcox3<- coxph(Surv(can.t_new, can.e)~  bpl.e_new  +
                  age +sex + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp,
                data=ds)

## Shizophrenia and overall cancer risk ----
ds<- ds %>% 
  mutate(sci.e_new = if_else(sci.t - can.t < -1, 1, 0),
         can.t_new = if_else(sci.e_new==1 & sci.t >0, can.t - sci.t, can.t))

ds$sci.e_new <- as.factor(ds$sci.e_new)

#model 1 
scicox1<- coxph(Surv(can.t_new, can.e)~  sci.e_new + 
                  age +sex + ethn + deprivation,
                data=ds)

#model 2
scicox2<- coxph(Surv(can.t_new, can.e)~  sci.e_new + 
                  age +sex + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                data=ds)

#model 3
scicox3<- coxph(Surv(can.t_new, can.e)~  sci.e_new + 
                  age +sex + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp ,
                data=ds)

## PTSD and Overall cancer risk -----

ds<- ds %>% 
  mutate(srd.pts.e_new = if_else(srd.pts.t - can.t < -1, 1, 0),
         can.t_new = if_else(srd.pts.e_new==1 & srd.pts.t >0, 
                             can.t - srd.pts.t, can.t))

ds$srd.pts.e_new <- as.factor(ds$srd.pts.e_new)

#model 1 
ptscox1<- coxph(Surv(can.t_new, can.e)~  srd.pts.e_new + 
                  age +sex + ethn + deprivation,
                data=ds)

#model 2
ptscox2<- coxph(Surv(can.t_new, can.e)~  srd.pts.e_new + 
                  age +sex + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                data=ds)

#model 3
ptscox3<- coxph(Surv(can.t_new, can.e)~  srd.pts.e_new + 
                  age +sex + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp ,
                data=ds)

## Testing the proportional hazard assumption for all models ----
phdep<- cox.zph(depcox3)
phanx<- cox.zph(anxcox3)
phbpl<- cox.zph(bplcox3)
phsci<- cox.zph(scicox3)
phpts<- cox.zph(ptscox3)

## Testing for interaction between Sex and MHCs -----
#Only models showing significant interaction were stratified by sex
ds <- uk %>% 
  mutate(dep.e_new = if_else(dep.t - can.t < -1 , 1, 0), 
         can.t_new = if_else(dep.e_new==1 & dep.t > 0, can.t - dep.t, can.t))

ds$dep.e_new<- as.factor(ds$dep.e_new)

interdep<- coxph(Surv(can.t_new, can.e) ~ 
                   dep.e_new + age +sex + dep.e_new*sex + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp , data=ds)

lrtest(interdep, depcox3)

# Depression model stratified by sex
fdepcox<- coxph(Surv(can.t_new, can.e) ~ dep.e_new +
                  age +sex + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp + menopause , 
                data= subset(ds, sex=="F")) 


mdepcox<- coxph(Surv(can.t_new, can.e) ~ dep.e_new +
                  age +sex + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp , data= subset(ds, sex=="M")) 

# Anxiety 
ds <- ds %>%
  mutate(anx.e_new = if_else(anx.t - can.t < -1, 1, 0),
         can.t_new = if_else(anx.e_new==1 & anx.t > 0, can.t - anx.t, can.t))

ds$anx.e_new<- as.factor(ds$anx.e_new)

interanx<- coxph(Surv(can.t_new, can.e) ~ 
                   anx.e_new + age +sex + anx.e_new*sex + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp , data=ds)

lrtest(interanx, anxcox3)

# Bipolar Disorders
ds<- ds %>%
  mutate(bpl.e_new = if_else(bpl.t - can.t < -1, 1, 0),
         can.t_new = if_else(bpl.e_new==1 & bpl.t >0, can.t - bpl.t, can.t))

ds$bpl.e_new<- factor(ds$bpl.e_new)

interbpl<- coxph(Surv(can.t_new, can.e) ~ 
                   bpl.e_new + age +sex + bpl.e_new*sex + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp , data=ds)

lrtest(interbpl, bplcox3)

# Bipolar Disorders model stratified by sex
fbplcox<- coxph(Surv(can.t_new, can.e)~  bpl.e_new +
                  age +sex + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp + menopause,
                data = subset(ds, sex=="F"))

mbplcox<- coxph(Surv(can.t_new, can.e)~  bpl.e_new +
                  age +sex + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp,
                data = subset(ds, sex=="M"))

# Schizophrenia 
ds<- ds %>%
  mutate(sci.e_new = if_else(sci.t - can.t < -1, 1, 0),
         can.t_new = if_else(sci.e_new==1 & sci.t >0, can.t - sci.t, can.t))

ds$sci.e_new <- as.factor(ds$sci.e_new)

intersci<- coxph(Surv(can.t_new, can.e)~  sci.e_new + 
                   age +sex + sci.e_new*sex + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp ,
                 data=ds)

lrtest(intersci, scicox3)

# PTSD

ds<- ds %>%
  mutate(srd.pts.e_new = if_else(srd.pts.t - can.t < -1, 1, 0),
         can.t_new = if_else(srd.pts.e_new==1 & srd.pts.t >0, 
                             can.t - srd.pts.t, can.t))

ds$srd.pts.e_new <- as.factor(ds$srd.pts.e_new)

interpts<- coxph(Surv(can.t_new, can.e)~  srd.pts.e_new + 
                   age +sex +srd.pts.e_new*sex+ ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp ,
                 data=ds)

lrtest(interpts, ptscox3)

## Cox Regression models with Site-specific cancer risk ----

### Depression and site-specific cancer risk ----

#Checking the number of incident cancer cases with depression to ensure powered analysis
crosstab(ds$can.bre.e, ds$dep.e, plot = F)
crosstab(ds$can.ova.e, ds$dep.e, plot = F)
crosstab(ds$can.ute.e, ds$dep.e, plot = F)
crosstab(ds$can.pro.e, ds$dep.e, plot = F)
crosstab(ds$can.lun.e, ds$dep.e, plot = F)
crosstab(ds$can.blo.e, ds$dep.e, plot = F)
crosstab(ds$can.col.e, ds$dep.e, plot = F)
crosstab(ds$can.liv.e, ds$dep.e, plot = F)

## Breast cancer 
ds <- ds %>% 
  filter(can.bre.t>0) %>%  
  mutate(dep.e_new = if_else(dep.t - can.bre.t < -3, 1, 0),
         can.bre.t_new = if_else(
           dep.e_new==1 & dep.t > 0, can.bre.t - dep.t, can.bre.t))

ds$dep.e_new<- factor(ds$dep.e_new)

#Model 1
dep.bre1<- coxph(Surv(can.bre.t_new, can.bre.e) ~ dep.e_new +
                   age+ ethn + deprivation , 
                 data= subset(ds, sex=="F")) 

#Model 2
dep.bre2<- coxph(Surv(can.bre.t_new, can.bre.e) ~ dep.e_new +
                   age  + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data= subset(ds, sex=="F")) 

#Model 3
dep.bre3<- coxph(Surv(can.bre.t_new, can.bre.e) ~ dep.e_new +
                   age  + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp + menopause , 
                 data= subset(ds, sex=="F")) 

#PH assumption
phdep.bre<- cox.zph(dep.bre3)

## Ovarian Cancer 
ds <- ds %>% 
  filter(can.ova.t>0) %>%  
  mutate(dep.e_new = if_else(dep.t - can.ova.t < -1, 1, 0),
         can.ova.t_new = if_else(
           dep.e_new==1 & dep.t > 0, can.ova.t - dep.t, can.ova.t))

ds$dep.e_new<- factor(ds$dep.e_new)

#Model 1
dep.ova1<- coxph(Surv(can.ova.t_new, can.ova.e) ~ dep.e_new + 
                   age +sex + ethn + deprivation , 
                 data= subset(ds, sex=="F")) 

#Model 2
dep.ova2<- coxph(Surv(can.ova.t_new, can.ova.e) ~ dep.e_new + 
                   age +sex + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data= subset(ds, sex=="F")) 

# Model 3
dep.ova3<- coxph(Surv(can.ova.t_new, can.ova.e) ~ dep.e_new + 
                   age +sex + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp + menopause , 
                 data= subset(ds, sex=="F")) 

#PH assumption
phdep.ova<- cox.zph(dep.ova3)

## Uterine Cancer 
ds <- ds %>% 
  filter(can.ute.t>0) %>%  
  mutate(dep.e_new = if_else(dep.t - can.ute.t < -1, 1, 0),
         can.ute.t_new = if_else(
           dep.e_new==1 & dep.t > 0, can.ute.t - dep.t, can.ute.t))

ds$dep.e_new<- factor(ds$dep.e_new)

#Model 1
dep.ute1<- coxph(Surv(can.ute.t_new, can.ute.e) ~dep.e_new +
                   age  + ethn + deprivation  , 
                 data= subset(ds, sex=="F")) 

#Model 2
dep.ute2<- coxph(Surv(can.ute.t_new, can.ute.e) ~dep.e_new +
                   age + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi , 
                 data= subset(ds, sex=="F")) 

#Model 3
dep.ute3<- coxph(Surv(can.ute.t_new, can.ute.e) ~dep.e_new +
                   age + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp + menopause , 
                 data= subset(ds, sex=="F")) 

#PH assumption 
phdep.ute<- cox.zph(dep.ute3)

## Lung Cancer 
ds <- ds %>% 
  filter(can.lun.t>0 ) %>%  
  mutate(dep.e_new = if_else(dep.t - can.lun.t < -1, 1, 0),
         can.lun.t_new = 
           if_else(dep.e_new==1 & dep.t > 0, can.lun.t - dep.t, can.lun.t))

ds$dep.e_new<- factor(ds$dep.e_new)

#Model 1
dep.lun1<- coxph(Surv(can.lun.t_new, can.lun.e) ~ dep.e_new + 
    age +sex + ethn + deprivation , data= ds)

#model 2
dep.lun2<- coxph(
  Surv(can.lun.t_new, can.lun.e) ~ dep.e_new + 
    age +sex + ethn + deprivation + 
    smok + alc.ut + met.tot + slp.durat + ffq.pmi, data= ds)

#model 3
dep.lun3<- coxph(
  Surv(can.lun.t_new, can.lun.e) ~ dep.e_new + 
    age +sex + ethn + deprivation + 
    smok + alc.ut + met.tot + slp.durat + ffq.pmi +
    bmi + blm.t2d + sbp, data= ds)

#PH Assumption
phdep.lun<- cox.zph(dep.lun3)

## Prostate Cancer
ds <- ds %>% 
  filter(can.pro.t>0) %>%  
  mutate(dep.e_new = if_else(dep.t - can.pro.t < -1, 1, 0),
         can.pro.t_new = 
           if_else(dep.e_new==1 & dep.t > 0, can.pro.t - dep.t, can.pro.t))

ds$dep.e_new<- factor(ds$dep.e_new)

#Model 1
dep.pro1<- coxph(Surv(can.pro.t_new, can.pro.e) ~ dep.e_new + 
                   age + ethn + deprivation ,
                 data= subset(ds, sex=="M"))

#Model 2
dep.pro2<- coxph(Surv(can.pro.t_new, can.pro.e) ~ dep.e_new + 
                   age + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data= subset(ds, sex=="M"))

#Model 3
dep.pro3<- coxph(Surv(can.pro.t_new, can.pro.e) ~ dep.e_new + 
                   age + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp,
                 data= subset(ds, sex=="M"))

#PH ASSUMPTION
phdep.pro<- cox.zph(dep.pro3)

## Blood cancer 
ds <- ds %>% 
  filter(can.blo.t>0) %>%  
  mutate(dep.e_new = if_else(dep.t - can.blo.t < -1, 1, 0),
         can.blo.t_new = 
           if_else(dep.e_new==1 & dep.t > 0, can.blo.t - dep.t, can.blo.t))

ds$dep.e_new<- factor(ds$dep.e_new)

#model 1
dep.blo1<- coxph(Surv(can.blo.t_new, can.blo.e) ~ dep.e_new +
                   age + sex+  ethn + deprivation , data= ds)

#model 2
dep.blo2<- coxph(Surv(can.blo.t_new, can.blo.e) ~ dep.e_new +
                   age + sex+  ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data= ds)

#model 3
dep.blo3<- coxph(Surv(can.blo.t_new, can.blo.e) ~ dep.e_new +
                   age + sex+  ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp , data= ds)

#PH ASSUMPTION
phdep.blo<- cox.zph(dep.blo3)

## Colon cancer 
ds <- ds %>% 
  filter(can.col.t>0) %>%  
  mutate(dep.e_new = if_else(dep.t - can.col.t < -1, 1, 0),
         can.col.t_new = 
           if_else(dep.e_new==1 & dep.t > 0, can.col.t - dep.t, can.col.t))

ds$dep.e_new<- factor(ds$dep.e_new)

#model 1
dep.col1<- coxph(Surv(can.col.t_new, can.col.e) ~ dep.e_new +
                   age + sex+  ethn + deprivation , data= ds)

#model 2
dep.col2<- coxph(Surv(can.col.t_new, can.col.e) ~ dep.e_new +
                   age + sex+  ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi
                 , data= ds)

#model 3
dep.col3<- coxph(Surv(can.col.t_new, can.col.e) ~ dep.e_new +
                   age + sex+  ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp , data= ds)

#PH Assumption
phdep.col<- cox.zph(dep.col3)

## Liver cancer 
ds <- ds %>% 
  filter(can.liv.t>0) %>%  
  mutate(dep.e_new = if_else(dep.t - can.liv.t < -1, 1, 0),
         can.liv.t_new = 
           if_else(dep.e_new==1 & dep.t > 0, can.liv.t - dep.t, can.liv.t))

ds$dep.e_new<- factor(ds$dep.e_new)

#model 1
dep.liv1<- coxph(Surv(can.liv.t_new, can.liv.e) ~ dep.e_new +
                   age + sex+  ethn + deprivation, data= ds)

#model 2
dep.liv2<- coxph(Surv(can.liv.t_new, can.liv.e) ~ dep.e_new +
                   age + sex+  ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data= ds)

#Model 3
dep.liv3<- coxph(Surv(can.liv.t_new, can.liv.e) ~ dep.e_new +
                   age + sex+  ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp, data= ds)

#PH assumption
phdep.liv<- cox.zph(dep.liv3)

### Anxiety and site-specific cancer risk ----
#Checking the number of incident cancer cases with anxiety  to ensure powered analysis
crosstab(ds$can.bre.e, ds$anx.e, plot = F)
crosstab(ds$can.ova.e, ds$anx.e, plot = F)
crosstab(ds$can.ute.e, ds$anx.e, plot = F)
crosstab(ds$can.pro.e, ds$anx.e, plot = F)
crosstab(ds$can.lun.e, ds$anx.e, plot = F)
crosstab(ds$can.blo.e, ds$anx.e, plot = F)
crosstab(ds$can.col.e, ds$anx.e, plot = F)
crosstab(ds$can.liv.e, ds$anx.e, plot = F)

## Breast cancer 
ds <- ds %>%  
  mutate(anx.e_new = if_else(anx.t - can.bre.t < -1, 1, 0),
         can.bre.t_new = if_else(
           anx.e_new==1 & anx.t > 0, can.bre.t - anx.t, can.bre.t))

ds$anx.e_new<- factor(ds$anx.e_new)

#model 1 
anx.bre1<- coxph(Surv(can.bre.t_new, can.bre.e) ~ anx.e_new + 
                   age  + ethn + deprivation,
                 data= subset(ds, sex=="F"))

#model 2
anx.bre2<- coxph(Surv(can.bre.t_new, can.bre.e) ~ anx.e_new + 
                   age  + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data= subset(ds, sex=="F"))

#model 3
anx.bre3<- coxph(Surv(can.bre.t_new, can.bre.e) ~ anx.e_new + 
                   age  + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp + menopause
                 , data= subset(ds, sex=="F"))

#PH Assumption
phanx.bre<- cox.zph(anx.bre3)

## Ovarian cancer 
ds <- ds %>%  
  mutate(anx.e_new = if_else(anx.t - can.bre.t < -1, 1, 0),
         can.ova.t_new = if_else(
           anx.e_new==1 & anx.t > 0, can.ova.t - anx.t, can.ova.t))

ds$anx.e_new<- factor(ds$anx.e_new)

#Model 1
anx.ova1<- coxph(Surv(can.ova.t_new, can.ova.e) ~ anx.e_new +
                   age  + ethn + deprivation,
                 data= subset(ds, sex=="F"))

#model 2
anx.ova2<- coxph(Surv(can.ova.t_new, can.ova.e) ~ anx.e_new +
                   age  + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data= subset(ds, sex=="F"))

#model 3
anx.ova3<- coxph(Surv(can.ova.t_new, can.ova.e) ~ anx.e_new +
                   age  + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp + menopause
                 , data= subset(ds, sex=="F"))

#PH Assumption
phanx.ova<- cox.zph(anx.ova3)

## Uterine Cancer 
ds <- ds %>%  
  mutate(anx.e_new = if_else(anx.t - can.ute.t < -1, 1, 0),
         can.ute.t_new = if_else(
           anx.e_new==1 & anx.t > 0, can.ute.t - anx.t, can.ute.t))

ds$anx.e_new<- factor(ds$anx.e_new)

#model 1
anx.ute1<- coxph(Surv(can.ute.t_new, can.ute.e) ~ anx.e_new +
                   age  + ethn + deprivation,
                 data= subset(ds, sex=="F"))

#model 2
anx.ute2<- coxph(Surv(can.ute.t_new, can.ute.e) ~ anx.e_new +
                   age  + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data= subset(ds, sex=="F"))

#model 3
anx.ute3<- coxph(Surv(can.ute.t_new, can.ute.e) ~ anx.e_new +
                   age  + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp + menopause,
                 data= subset(ds, sex=="F"))

#PH Assumption
phanx.ute<- cox.zph(anx.ute3)

## Lung cancer 
ds <- ds %>%  
  mutate(anx.e_new = if_else(anx.t - can.lun.t < -1, 1, 0),
         can.lun.t_new = if_else(
           anx.e_new==1 & anx.t > 0, can.lun.t - anx.t, can.lun.t))

ds$anx.e_new<- factor(ds$anx.e_new)

#Model 1
anx.lun1<- coxph(Surv(can.lun.t_new, can.lun.e) ~ anx.e_new + 
                   age  + sex + ethn + deprivation,
                 data= ds)

#model 2
anx.lun2<- coxph(Surv(can.lun.t_new, can.lun.e) ~ anx.e_new + 
                   age  + sex + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data= ds)

#Model 3
anx.lun3<- coxph(Surv(can.lun.t_new, can.lun.e) ~ anx.e_new + 
                   age  + sex + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp,
                 data= ds)

#PH assumption
phanx.lun<- cox.zph(anx.lun3)

## Prostate Cancer 
ds <- ds %>% 
  mutate(anx.e_new = if_else(anx.t - can.pro.t < -1, 1, 0),
         can.pro.t_new = if_else(
           anx.e_new==1 & anx.t > 0, can.pro.t - anx.t, can.pro.t))

ds$anx.e_new<- factor(ds$anx.e_new)

#Model 1
anx.pro1<- coxph(Surv(can.pro.t_new, can.pro.e) ~ anx.e_new + 
                   age + ethn + deprivation,
                 data= subset(ds, sex=="M"))

#Model 2
anx.pro2<- coxph(Surv(can.pro.t_new, can.pro.e) ~ anx.e_new + 
                   age + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data= subset(ds, sex=="M"))

#Model 3
anx.pro3<- coxph(Surv(can.pro.t_new, can.pro.e) ~ anx.e_new + 
                   age + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp,
                 data= subset(ds, sex=="M"))

#PH assumption
phanx.pro<- cox.zph(anx.pro3)


## Blood Cancer
ds <- ds %>%  
  mutate(anx.e_new = if_else(anx.t - can.blo.t < -1, 1, 0),
         can.blo.t_new = if_else(
           anx.e_new==1 & anx.t > 0, can.blo.t - anx.t, can.blo.t))

ds$anx.e_new<- factor(ds$anx.e_new)

#Model 1
anx.blo1<- coxph(Surv(can.blo.t_new, can.blo.e) ~ anx.e_new + 
                   age + sex+ ethn + deprivation, data= ds)

#Model 2
anx.blo2<- coxph(Surv(can.blo.t_new, can.blo.e) ~ anx.e_new + 
                   age + sex+ ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi, data= ds)

#Model 3
anx.blo3<- coxph(Surv(can.blo.t_new, can.blo.e) ~ anx.e_new + 
                   age + sex+ ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp, data= ds)

#Ph assumption
phanx.blo<- cox.zph(anx.blo3)

## Colon cancer 
ds<- ds %>%
  mutate(anx.e_new = if_else(anx.t - can.blo.t < -1, 1, 0),
         can.col.t_new = if_else(
           anx.e_new==1 & anx.t > 0, can.col.t - anx.t, can.col.t))

ds$anx.e_new<- factor(ds$anx.e_new)

#model 1
anx.col1<- coxph(Surv(can.col.t_new, can.col.e) ~ anx.e_new + 
                   age + sex+ ethn + deprivation , 
                 data= ds)

#model 2
anx.col2<- coxph(Surv(can.col.t_new, can.col.e) ~ anx.e_new + 
                   age + sex+ ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi, 
                 data= ds)

#Model 3
anx.col3<- coxph(Surv(can.col.t_new, can.col.e) ~ anx.e_new + 
                   age + sex+ ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp, 
                 data= ds)

#PH assumption
phanx.col<- cox.zph(anx.col3)

## Liver Cancer 
ds<- ds %>% 
  filter(can.liv.t>0) %>%  
  mutate(anx.e_new = if_else(anx.t - can.liv.t < -1, 1, 0),
         can.liv.t_new = if_else(
           anx.e_new==1 & anx.t > 0, can.liv.t - anx.t, can.liv.t))

ds$anx.e_new<- factor(ds$anx.e_new)

#Model 1
anx.liv1<- coxph(Surv(can.liv.t_new, can.liv.e) ~ anx.e_new + 
                   age + sex+ ethn + deprivation, 
                 data= ds)

#Model 2
anx.liv2<- coxph(Surv(can.liv.t_new, can.liv.e) ~ anx.e_new + 
                   age + sex+ ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi , 
                 data= ds)

#Model 3
anx.liv3<- coxph(Surv(can.liv.t_new, can.liv.e) ~ anx.e_new + 
                   age + sex+ ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp, 
                 data= ds)

summary(anx.liv)

#PH assumption
phanx.liv<- cox.zph(anx.liv3)

### Bipolar Disorder and site-specific cancer risk ----
#Checking the number of incident cancer cases with BD to ensure powered analysis

crosstab(ds$can.bre.e, ds$bpl.e, plot = F)
crosstab(ds$can.ova.e, ds$bpl.e, plot = F)
crosstab(ds$can.ute.e, ds$bpl.e, plot = F)
crosstab(ds$can.pro.e, ds$bpl.e, plot = F)
crosstab(ds$can.lun.e, ds$bpl.e, plot = F)
crosstab(ds$can.blo.e, ds$bpl.e, plot = F)
crosstab(ds$can.col.e, ds$bpl.e, plot = F)
crosstab(ds$can.liv.e, ds$bpl.e, plot = F)

## Breast cancer 
ds <- ds %>% 
  mutate(bpl.e_new = if_else(bpl.t - can.bre.t < -1, 1, 0),
         can.bre.t_new = if_else(bpl.e_new==1 & bpl.t >0,
                                 can.bre.t - bpl.t, can.bre.t))

ds$bpl.e_new<- factor(ds$bpl.e_new)

#model 1
bpl.bre1<- coxph(Surv(can.bre.t_new, can.bre.e)~  bpl.e_new +
                   age + ethn + deprivation, 
                 data= subset(ds, sex=="F"))

#model 2
bpl.bre2<- coxph(Surv(can.bre.t_new, can.bre.e)~  bpl.e_new +
                   age + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi, 
                 data= subset(ds, sex=="F"))

#Model 3
bpl.bre3<- coxph(Surv(can.bre.t_new, can.bre.e)~  bpl.e_new +
                   age + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp + menopause, 
                 data= subset(ds, sex=="F"))

#PH assumption 
phbpl.bre<- cox.zph(bpl.bre3)

## Prostate cancer 
ds <- ds %>% 
  mutate(bpl.e_new = if_else(bpl.t - can.pro.t < -1, 1, 0),
         can.pro.t_new = if_else(bpl.e_new==1 & bpl.t >0,
                                 can.pro.t - bpl.t, can.pro.t))

ds$bpl.e_new<- factor(ds$bpl.e_new)

#model 1
bpl.pro1<- coxph(Surv(can.pro.t_new, can.pro.e)~  bpl.e_new  +
                   age + ethn + deprivation , 
                 data=subset(ds, sex== "M"))

#model 2
bpl.pro2<- coxph(Surv(can.pro.t_new, can.pro.e)~  bpl.e_new  +
                   age + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi, 
                 data=subset(ds, sex== "M"))

#model 3
bpl.pro3<- coxph(Surv(can.pro.t_new, can.pro.e)~  bpl.e_new  +
                   age + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp , 
                 data=subset(ds, sex== "M"))

#PH assumption
phbpl.pro<- cox.zph(bpl.pro3)

## Blood cancer
ds<- ds %>% 
  mutate(bpl.e_new = if_else(bpl.t - can.blo.t < -1, 1, 0),
         can.blo.t_new = if_else(bpl.e_new==1 & bpl.t >0,
                                 can.blo.t - bpl.t, can.blo.t))

ds$bpl.e_new<- factor(ds$bpl.e_new)

#model 1
bpl.blo1<- coxph(Surv(can.blo.t_new, can.blo.e)~  bpl.e_new +
                   age + sex + ethn + deprivation, 
                 data= ds)

#model 2
bpl.blo2<- coxph(Surv(can.blo.t_new, can.blo.e)~  bpl.e_new +
                   age + sex + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi, 
                 data= ds)

#model 3
bpl.blo3<- coxph(Surv(can.blo.t_new, can.blo.e)~  bpl.e_new +
                   age + sex + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp, 
                 data= ds)

#PH ASSUMPTION 
phbpl.blo<- cox.zph(bpl.blo3)



## Sensitivity analysis -----

#Standardisation of PHQ-4 scores 
ds$phqs<- scale(ds$phq)

##Running models with PHQ-4 scores as the main independent variable 

#With overall cancer risk
sensall<- coxph(Surv(can.t, can.e) ~  phqs+
                  age +sex + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp , data=ds)

## With site-specific cancer risk
# Breast cancer
sensbre<- coxph(Surv(can.bre.t, can.bre.e) ~ 
                  phqs +age + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp + menopause 
                , data= subset(ds, sex=="F"))

#ovarian 
sensova<- coxph(Surv(can.ova.t, can.ova.e) ~ 
                  phqs +age + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp + menopause 
                , data= subset(uk, sex=="F"))

#uterine 
sensute<- coxph(Surv(can.ute.t, can.ute.e) ~ 
                  phqs +age + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp + menopause 
                , data= subset(uk, sex=="F"))

#Prostate cancer 
senspro<- coxph(Surv(can.pro.t, can.pro.e) ~ 
                  phqs +age + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp 
                , data= subset(uk, sex=="M"))

# lung cancer 
senslun<- coxph(Surv(can.lun.t, can.lun.e) ~ phqs+ 
                  age +sex + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp , data=ds)

#Blood Cancer 
sensblo<- coxph(Surv(can.blo.t, can.blo.e) ~ 
                  phqs +age + sex+ ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp , data= uk)

#Colon Cancer 
senscol<- coxph(Surv(can.col.t, can.col.e) ~ 
                  phqs +age + sex+  ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp , data= uk)

#Liver Cancer 
sensliv<- coxph(Surv(can.liv.t, can.liv.e) ~ 
                  phqs +age + sex+  ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp , data= uk)
