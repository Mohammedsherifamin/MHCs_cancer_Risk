#rm(list = ls())
source("scripts/1.data_prepration.R")

#Loading relevant packages
library(tidyverse)
library(descr)
library(mosaic)
library(broom)
library(Hmisc)
library(lmtest)
library(survival)
library(gtExtras)

#select only relevant variables for the main analysis
uk <- uk %>%
  dplyr::select(eid, age, sex,ethn,
                deprivation, smok,alc.ut, slp.durat, 
                met.tot, ffq.pmi, bmi, blm.t2d, sbp, menopause,
                dep.e, anx.e, bpl.e, sci.e, srd.pts.e, dep.t, anx.t, bpl.t, 
                sci.t, srd.pts.t, can.e, can.t, can.bre.e, can.ova.e, can.ute.e, 
                can.pro.e, can.lun.e, can.blo.e, can.col.e, can.liv.e, can.bre.t, 
                can.ova.t, can.ute.t, can.pro.t, can.lun.t ,can.blo.t , 
                can.col.t, can.liv.t)

nrow(uk)

## Due to the small number of missing data, we decided to drop missing data and preform a 
#complete cases analysis
uk <- uk %>%  drop_na()

nrow(uk)

# Cox Regression Models with overall cancer ------
## Depression and overall cancer risk-----

#Create a new data set
ds <- uk %>% 
  # a new variable to ensure a 1 year land mark analysis and that post baseline depression cases
  # are only those with depression diagnosis before cancer diagnosis 
  mutate(dep.e_new = if_else(dep.t - can.t < -1 , 1, 0), 
         can.t_new = if_else(dep.e_new==1 & dep.t > 0, can.t - dep.t, can.t))

ds$dep.e_new<- as.factor(ds$dep.e_new)

#Adjusted for Demographics "model 1"
depcox1<- coxph(Surv(can.t_new, can.e)~ dep.e_new +
                  age +sex + ethn + deprivation, 
                data = ds)

#Adjusted for Sociodemographics + health conditions "model 2"
depcox2<- coxph(Surv(can.t_new, can.e)~ dep.e_new +
                  age + sex + ethn + deprivation + 
                  bmi + blm.t2d + sbp,
                data = ds)

#Fully Adjusted "model 3"
depcox3 <- coxph(Surv(can.t_new, can.e) ~ dep.e_new + 
                   age + sex + ethn + deprivation + 
                   bmi + blm.t2d + sbp +
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data = ds)

#Testing for PH assumption
phdep<- cox.zph(depcox3)

## Anxiety and overall cancer risk ----
ds <- uk %>% 
  mutate(anx.e_new = if_else(anx.t - can.t < -1, 1, 0),
         can.t_new = if_else(anx.e_new==1 & anx.t > 0, can.t - anx.t, can.t))

ds$anx.e_new<- as.factor(ds$anx.e_new)

#model 1
anxcox1<- coxph(Surv(can.t_new, can.e) ~ anx.e_new+ 
                  age +sex + ethn + deprivation, 
                data=ds) 

#model 2 
anxcox2<- coxph(Surv(can.t_new, can.e) ~ anx.e_new + 
                  age + sex + ethn + deprivation + 
                  bmi + blm.t2d + sbp,
                data = ds)

#model 3
anxcox3<- coxph(Surv(can.t_new, can.e) ~ anx.e_new + 
                  age + sex + ethn + deprivation + 
                  bmi + blm.t2d + sbp +
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                data = ds)

#Testing for PH assumption
phanx<- cox.zph(anxcox3)


## Bipolar Disorder and overall cancer risk -----
ds<- uk %>%  
  mutate(bpl.e_new = if_else(bpl.t - can.t < -1, 1, 0),
         can.t_new = if_else(bpl.e_new==1 & bpl.t >0, can.t - bpl.t, can.t))

ds$bpl.e_new<- factor(ds$bpl.e_new)

#model 1
bplcox1<- coxph(Surv(can.t_new, can.e)~  bpl.e_new  + age +sex + ethn + 
                  deprivation , data=ds)

#model 2
bplcox2<- coxph(Surv(can.t_new, can.e)~  bpl.e_new  +
                  age + sex + ethn + deprivation + 
                  bmi + blm.t2d + sbp,
                data = ds)
#model 3
bplcox3<- coxph(Surv(can.t_new, can.e)~  bpl.e_new  +
                  age + sex + ethn + deprivation + 
                  bmi + blm.t2d + sbp +
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                data = ds)

#Testing for PH assumption
phbpl<- cox.zph(bplcox3)


## Shizophrenia and overall cancer risk ----
ds<- uk %>% 
  mutate(sci.e_new = if_else(sci.t - can.t < -1, 1, 0),
         can.t_new = if_else(sci.e_new==1 & sci.t >0, can.t - sci.t, can.t))

ds$sci.e_new <- as.factor(ds$sci.e_new)

#model 1 
scicox1<- coxph(Surv(can.t_new, can.e)~  sci.e_new + 
                  age +sex + ethn + deprivation,
                data=ds)

#model 2
scicox2<- coxph(Surv(can.t_new, can.e)~  sci.e_new + 
                  age + sex + ethn + deprivation + 
                  bmi + blm.t2d + sbp,
                data = ds)

#model 3
scicox3<- coxph(Surv(can.t_new, can.e)~  sci.e_new + 
                  age + sex + ethn + deprivation + 
                  bmi + blm.t2d + sbp +
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                data = ds)

#Testing for PH assumption
phsci<- cox.zph(scicox3)

## PTSD and Overall cancer risk -----

ds<- uk %>% 
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
                  age + sex + ethn + deprivation + 
                  bmi + blm.t2d + sbp,
                data = ds)

#model 3
ptscox3<- coxph(Surv(can.t_new, can.e)~  srd.pts.e_new + 
                  age + sex + ethn + deprivation + 
                  bmi + blm.t2d + sbp +
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                data = ds)

# Testing the PH assumption
phpts<- cox.zph(ptscox1)





# Interaction between Sex and MHCs #################################################

#Only models showing significant interaction were stratified by sex
ds <- uk %>% 
  mutate(dep.e_new = if_else(dep.t - can.t < -1 , 1, 0), 
         can.t_new = if_else(dep.e_new==1 & dep.t > 0, can.t - dep.t, can.t))

ds$dep.e_new<- as.factor(ds$dep.e_new)

interdep<- coxph(Surv(can.t_new, can.e) ~ 
                   dep.e_new + age +sex + dep.e_new*sex + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp , data=ds)

int_dep<- lrtest(interdep, depcox3)

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
ds <- uk %>%
  mutate(anx.e_new = if_else(anx.t - can.t < -1, 1, 0),
         can.t_new = if_else(anx.e_new==1 & anx.t > 0, can.t - anx.t, can.t))

ds$anx.e_new<- as.factor(ds$anx.e_new)

interanx<- coxph(Surv(can.t_new, can.e) ~ 
                   anx.e_new + age +sex + anx.e_new*sex + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp , data=ds)

int_anx<- lrtest(interanx, anxcox3)

# Bipolar Disorders
ds<- uk %>%
  mutate(bpl.e_new = if_else(bpl.t - can.t < -1, 1, 0),
         can.t_new = if_else(bpl.e_new==1 & bpl.t >0, can.t - bpl.t, can.t))

ds$bpl.e_new<- factor(ds$bpl.e_new)

interbpl<- coxph(Surv(can.t_new, can.e) ~ 
                   bpl.e_new + age +sex + bpl.e_new*sex + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp , data=ds)

int_bpl <- lrtest(interbpl, bplcox3)

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
ds<- uk %>%
  mutate(sci.e_new = if_else(sci.t - can.t < -1, 1, 0),
         can.t_new = if_else(sci.e_new==1 & sci.t >0, can.t - sci.t, can.t))

ds$sci.e_new <- as.factor(ds$sci.e_new)

intersci<- coxph(Surv(can.t_new, can.e)~  sci.e_new + 
                   age +sex + sci.e_new*sex + ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp ,
                 data=ds)

int_sci<- lrtest(intersci, scicox3)

# PTSD
ds<- uk %>%
  mutate(srd.pts.e_new = if_else(srd.pts.t - can.t < -1, 1, 0),
         can.t_new = if_else(srd.pts.e_new==1 & srd.pts.t >0, 
                             can.t - srd.pts.t, can.t))

ds$srd.pts.e_new <- as.factor(ds$srd.pts.e_new)

interpts<- coxph(Surv(can.t_new, can.e)~  srd.pts.e_new + 
                   age +sex +srd.pts.e_new*sex+ ethn + deprivation + 
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                   bmi + blm.t2d + sbp ,
                 data=ds)

int_pts<- lrtest(interpts, ptscox3)


## Table for interaction results -----

library(broom)
library(dplyr)
library(gt)

# Create a named list of your interaction test results
lrt_list <- list(
  "Depression" = int_dep,
  "Anxiety" = int_anx,
  "Bipolar disorder" = int_bpl,
  "SZ" = int_sci,
  "PTS" = int_pts
)

# Extract the 2nd row (comparison between full and reduced model)
interaction_results <- purrr::map_dfr(
  .x = lrt_list,
  .f = ~ tidy(.x)[2, ],
  .id = "Interaction Term") %>%
  select(`Interaction Term`, statistic, df, p.value) %>%
  rename(LRT = statistic, `df` = df, `P-value` = p.value) %>%
  mutate(
    LRT = round(LRT, 3),
    `P-value` = ifelse(`P-value` < 0.05,
                       paste0(formatC(`P-value`, digits = 4, format = "f"), "*"),
                       formatC(`P-value`, digits = 4, format = "f")))


tbl.int <- flextable(interaction_results) %>%
  set_header_labels(
    `Interaction Term` = "Interaction Term",
    LRT = "LRT¹",
    df = "df",
    `P-value` = "P-value"
  ) %>%
  add_footer_lines("¹Likelihood Ratio Tests comparing models with and without a sex interaction term for each condition. All models are adjusted for age, ethnicity, deprivation, smoking status, alcohol intake, physical activity, sleep duration, processed meat intake, BMI, diabetes, and systolic blood pressure.") %>%
  autofit()

tbl.int %>% save_as_docx(path = "output/tbl.int.docx")



# Site-specific cancer risk #################################################

#Checking the number of incident cancer cases with MHCs to ensure powered analysis
#  Table for MHCs and site-specific cancer ---------
#Adding the overall MHCs variable again 

uk <-uk %>% mutate(mhc.e= 
                     if_else(dep.e==1 | anx.e ==1 | bpl.e==1 |
                               srd.pts.e == 1 | sci.e== 1, 1,0))

#Table for overall
tbl.mhc <- tbl_summary(
  uk %>%
    filter(mhc.e == 1) %>%
    mutate(mhc.e = case_when(mhc.e == 1 ~ "All MHCs")) %>%
    select(
      mhc.e, can.e,
      can.bre.e, can.ova.e, can.ute.e, can.pro.e,
      can.lun.e, can.blo.e, can.col.e, can.liv.e),
  by = mhc.e,
  label = list(
    can.e    ~ "Incidence of Overall",
    can.bre.e ~ "Breast",
    can.ova.e ~ "Ovarian",
    can.ute.e ~ "Uterine",
    can.pro.e ~ "Prostate",
    can.lun.e ~ "Lung",
    can.blo.e ~ "Blood",
    can.col.e ~ "Colorectal",
    can.liv.e ~ "Liver")) %>%
  modify_header(label = "**Cancer incidence**") 


#Table for depression
tbl.dep <- tbl_summary(
  uk %>%
    filter(dep.e == 1) %>%
    mutate(dep.e = case_when(dep.e == 1 ~ "DD")) %>%
    select(
      dep.e, can.e,
      can.bre.e, can.ova.e, can.ute.e, can.pro.e,
      can.lun.e, can.blo.e, can.col.e, can.liv.e),
  by = dep.e,
  label = list(
    can.e    ~ "Incidence of Overall",
    can.bre.e ~ "Breast",
    can.ova.e ~ "Ovarian",
    can.ute.e ~ "Uterine",
    can.pro.e ~ "Prostate",
    can.lun.e ~ "Lung",
    can.blo.e ~ "Blood",
    can.col.e ~ "Colorectal",
    can.liv.e ~ "Liver")) %>%
  modify_header(label = "**Cancer incidence**") 


#Table for anxiety disorders
tbl.anx <- tbl_summary(
  uk %>%
    filter(anx.e == 1) %>%
    mutate(anx.e = case_when(anx.e == 1 ~ "AD")) %>%
    select(
      anx.e, can.e, 
      can.bre.e, can.ova.e, can.ute.e, can.pro.e,
      can.lun.e, can.blo.e, can.col.e, can.liv.e),
  by = anx.e,
  label = list(
    can.e    ~ "Incidence of Overall",
    can.bre.e ~ "Breast",
    can.ova.e ~ "Ovarian",
    can.ute.e ~ "Uterine",
    can.pro.e ~ "Prostate",
    can.lun.e ~ "Lung",
    can.blo.e ~ "Blood",
    can.col.e ~ "Colorectal",
    can.liv.e ~ "Liver")) %>%
  modify_header(label = "**Cancer incidence**") 

#Table for bipolar disorders
tbl.bpl <- tbl_summary(
  uk %>%
    filter(bpl.e == 1) %>%
    mutate(bpl.e = case_when(bpl.e == 1 ~ "BD")) %>%
    select(
      bpl.e, can.e,
      can.bre.e, can.ova.e, can.ute.e, can.pro.e,
      can.lun.e, can.blo.e, can.col.e, can.liv.e),
  by = bpl.e,
  label = list(
    can.e    ~ "Incidence of Overall",
    can.bre.e ~ "Breast",
    can.ova.e ~ "Ovarian",
    can.ute.e ~ "Uterine",
    can.pro.e ~ "Prostate",
    can.lun.e ~ "Lung",
    can.blo.e ~ "Blood",
    can.col.e ~ "Colorectal",
    can.liv.e ~ "Liver")) %>%
  modify_header(label = "**Cancer incidence**") 

#table for SZ
tbl.sci <- tbl_summary(
  uk %>%
    filter(sci.e == 1) %>%
    mutate(sci.e = case_when(sci.e == 1 ~ "SZ")) %>%
    select(
      sci.e, can.e,
      can.bre.e, can.ova.e, can.ute.e, can.pro.e,
      can.lun.e, can.blo.e, can.col.e, can.liv.e),
  by = sci.e,
  label = list(
    can.e    ~ "Incidence of Overall",
    can.bre.e ~ "Breast",
    can.ova.e ~ "Ovarian",
    can.ute.e ~ "Uterine",
    can.pro.e ~ "Prostate",
    can.lun.e ~ "Lung",
    can.blo.e ~ "Blood",
    can.col.e ~ "Colorectal",
    can.liv.e ~ "Liver")) %>%
  modify_header(label = "**Cancer incidence**") 

#Table for ptsd
tbl.pts <- tbl_summary(
  uk %>%
    filter(srd.pts.e == 1) %>%
    mutate(srd.pts.e = case_when(srd.pts.e == 1 ~ "PTSD")) %>%
    select(
      srd.pts.e, can.e,
      can.bre.e, can.ova.e, can.ute.e, can.pro.e,
      can.lun.e, can.blo.e, can.col.e, can.liv.e),
  by = srd.pts.e,
  label = list(
    can.e    ~ "Incidence of Overall",
    can.bre.e ~ "Breast",
    can.ova.e ~ "Ovarian",
    can.ute.e ~ "Uterine",
    can.pro.e ~ "Prostate",
    can.lun.e ~ "Lung",
    can.blo.e ~ "Blood",
    can.col.e ~ "Colorectal",
    can.liv.e ~ "Liver")) %>%
  modify_header(label = "**Cancer incidence**") 


#Merging all tables and saving the output
tbl_merge(tbls = list(tbl.mhc, tbl.dep, tbl.anx, tbl.bpl, tbl.sci, tbl.pts)) %>%
  as_flex_table() %>%
  add_footer_lines(
    "MHCs: mental health conditions; DD: depressive disorders; AD: anxiety disorders; BD: bipolar disorders; SZ: schizophrenia; PTSD: post-traumatic stress disorders"
  ) %>%
  save_as_docx(path = "output/tbl_cancer_incidence.docx")


#Only DD and AD meet the requirement of 5 cases per predictor variables.


## Depression -----
### Breast cancer -----
ds <- uk %>% 
  filter(can.bre.t>0) %>%  
  mutate(dep.e_new = if_else(dep.t - can.bre.t < -1, 1, 0),
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
                   bmi + blm.t2d + sbp + menopause, 
                 data= subset(ds, sex=="F")) 

#Model 3
dep.bre3<- coxph(Surv(can.bre.t_new, can.bre.e) ~ dep.e_new +
                   age  + ethn + deprivation + 
                   bmi + blm.t2d + sbp + menopause +
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data= subset(ds, sex=="F")) 

#PH assumption
phdep.bre<- cox.zph(dep.bre3)

### Ovarian Cancer ------
ds <- uk %>% 
  filter(can.ova.t>0) %>%  
  mutate(dep.e_new = if_else(dep.t - can.ova.t < -1, 1, 0),
         can.ova.t_new = if_else(
           dep.e_new==1 & dep.t > 0, can.ova.t - dep.t, can.ova.t))

ds$dep.e_new<- factor(ds$dep.e_new)

#Model 1
dep.ova1<- coxph(Surv(can.ova.t_new, can.ova.e) ~ dep.e_new + 
                   age + ethn + deprivation , 
                 data= subset(ds, sex=="F")) 

#Model 2
dep.ova2<- coxph(Surv(can.ova.t_new, can.ova.e) ~ dep.e_new + 
                   age  + ethn + deprivation + 
                   bmi + blm.t2d + sbp + menopause, 
                 data= subset(ds, sex=="F")) 

# Model 3
dep.ova3<- coxph(Surv(can.ova.t_new, can.ova.e) ~ dep.e_new + 
                   age  + ethn + deprivation + 
                   bmi + blm.t2d + sbp + menopause +
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data= subset(ds, sex=="F")) 

#PH assumption
phdep.ova<- cox.zph(dep.ova3)

### Uterine Cancer -----
ds <- uk %>% 
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
                   age  + ethn + deprivation + 
                   bmi + blm.t2d + sbp + menopause, 
                 data= subset(ds, sex=="F")) 

#Model 3
dep.ute3<- coxph(Surv(can.ute.t_new, can.ute.e) ~dep.e_new +
                   age  + ethn + deprivation + 
                   bmi + blm.t2d + sbp + menopause +
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data= subset(ds, sex=="F")) 

#PH assumption 
phdep.ute<- cox.zph(dep.ute3)

### Prostate Cancer-----
ds <- uk %>% 
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
                   bmi + blm.t2d + sbp,
                 data= subset(ds, sex=="M"))

#Model 3
dep.pro3<- coxph(Surv(can.pro.t_new, can.pro.e) ~ dep.e_new + 
                   age +  ethn + deprivation + 
                   bmi + blm.t2d + sbp +
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data= subset(ds, sex=="M"))

#PH ASSUMPTION
phdep.pro<- cox.zph(dep.pro3)

### Lung Cancer ------
ds <- uk %>% 
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
    age + sex + ethn + deprivation + 
    bmi + blm.t2d + sbp,
  data = ds)

#model 3
dep.lun3<- coxph(
  Surv(can.lun.t_new, can.lun.e) ~ dep.e_new + 
    age + sex + ethn + deprivation + 
    bmi + blm.t2d + sbp +
    smok + alc.ut + met.tot + slp.durat + ffq.pmi,
  data = ds)

#PH Assumption
phdep.lun<- cox.zph(dep.lun3)


### Blood cancer ------
ds <- uk %>% 
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
                   age + sex + ethn + deprivation + 
                   bmi + blm.t2d + sbp,
                 data = ds)

#model 3
dep.blo3<- coxph(Surv(can.blo.t_new, can.blo.e) ~ dep.e_new +
                   age + sex + ethn + deprivation + 
                   bmi + blm.t2d + sbp +
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data = ds)

#PH ASSUMPTION
phdep.blo<- cox.zph(dep.blo3)

### Colon cancer ------
ds <- uk %>% 
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
                   age + sex + ethn + deprivation + 
                   bmi + blm.t2d + sbp,
                 data = ds)

#model 3
dep.col3<- coxph(Surv(can.col.t_new, can.col.e) ~ dep.e_new +
                   age + sex + ethn + deprivation + 
                   bmi + blm.t2d + sbp +
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data = ds)
#PH Assumption
phdep.col<- cox.zph(dep.col3)

### Liver cancer ------
ds <- uk %>% 
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
                   age + sex + ethn + deprivation + 
                   bmi + blm.t2d + sbp,
                 data = ds)

#Model 3
dep.liv3<- coxph(Surv(can.liv.t_new, can.liv.e) ~ dep.e_new +
                   age + sex + ethn + deprivation + 
                   bmi + blm.t2d + sbp +
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data = ds)

#PH assumption
phdep.liv<- cox.zph(dep.liv3)


## Anxiety ------

### Breast cancer ----

ds <- uk %>%  
  filter(can.bre.t>0) %>% 
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
                   bmi + blm.t2d + sbp + menopause, 
                 data= subset(ds, sex=="F")) 

#model 3
anx.bre3<- coxph(Surv(can.bre.t_new, can.bre.e) ~ anx.e_new + 
                   age  + ethn + deprivation + 
                   bmi + blm.t2d + sbp + menopause +
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data= subset(ds, sex=="F")) 

#PH Assumption
phanx.bre<- cox.zph(anx.bre3)


### Ovarian cancer -----
ds <- uk %>% 
  filter(can.ova.t>0) %>%  
  mutate(anx.e_new = if_else(anx.t - can.ova.t < -1, 1, 0),
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
                   bmi + blm.t2d + sbp + menopause, 
                 data= subset(ds, sex=="F")) 

#model 3
anx.ova3<- coxph(Surv(can.ova.t_new, can.ova.e) ~ anx.e_new +
                   age  + ethn + deprivation + 
                   bmi + blm.t2d + sbp + menopause +
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data= subset(ds, sex=="F")) 

#PH Assumption
phanx.ova<- cox.zph(anx.ova3)

### Uterine Cancer ------
ds <- uk %>% 
  filter(can.ute.t>0) %>%  
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
                   bmi + blm.t2d + sbp + menopause, 
                 data= subset(ds, sex=="F")) 

#model 3
anx.ute3<- coxph(Surv(can.ute.t_new, can.ute.e) ~ anx.e_new +
                   age  + ethn + deprivation + 
                   bmi + blm.t2d + sbp + menopause +
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data= subset(ds, sex=="F")) 

#PH Assumption
phanx.ute<- cox.zph(anx.ute3)

### Prostate Cancer -----
ds <- uk %>% 
  filter(can.pro.t>0) %>% 
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
                   bmi + blm.t2d + sbp,
                 data= subset(ds, sex=="M"))

#Model 3
anx.pro3<- coxph(Surv(can.pro.t_new, can.pro.e) ~ anx.e_new + 
                   age +  ethn + deprivation + 
                   bmi + blm.t2d + sbp +
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data= subset(ds, sex=="M"))

#PH assumption
phanx.pro<- cox.zph(anx.pro3)


### Lung cancer -------
ds <- uk %>%  
  filter(can.lun.t>0) %>% 
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
                   age + sex + ethn + deprivation + 
                   bmi + blm.t2d + sbp,
                 data = ds)

#Model 3
anx.lun3<- coxph(Surv(can.lun.t_new, can.lun.e) ~ anx.e_new + 
                   age + sex + ethn + deprivation + 
                   bmi + blm.t2d + sbp +
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data = ds)

#PH assumption
phanx.lun<- cox.zph(anx.lun3)


### Blood Cancer--------
ds <- uk %>% 
  filter(can.blo.t>0) %>%  
  mutate(anx.e_new = if_else(anx.t - can.blo.t < -1, 1, 0),
         can.blo.t_new = if_else(
           anx.e_new==1 & anx.t > 0, can.blo.t - anx.t, can.blo.t))

ds$anx.e_new<- factor(ds$anx.e_new)

#Model 1
anx.blo1<- coxph(Surv(can.blo.t_new, can.blo.e) ~ anx.e_new + 
                   age + sex+ ethn + deprivation, data= ds)

#Model 2
anx.blo2<- coxph(Surv(can.blo.t_new, can.blo.e) ~ anx.e_new + 
                   age + sex + ethn + deprivation + 
                   bmi + blm.t2d + sbp,
                 data = ds)

#Model 3
anx.blo3<- coxph(Surv(can.blo.t_new, can.blo.e) ~ anx.e_new + 
                   age + sex + ethn + deprivation + 
                   bmi + blm.t2d + sbp +
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data = ds)

#Ph assumption
phanx.blo<- cox.zph(anx.blo3)

### Colon cancer ----
ds<- uk %>%
  filter(can.col.t >0) %>% 
  mutate(anx.e_new = if_else(anx.t - can.col.t < -1, 1, 0),
         can.col.t_new = if_else(
           anx.e_new==1 & anx.t > 0, can.col.t - anx.t, can.col.t))

ds$anx.e_new<- factor(ds$anx.e_new)

#model 1
anx.col1<- coxph(Surv(can.col.t_new, can.col.e) ~ anx.e_new + 
                   age + sex+ ethn + deprivation , 
                 data= ds)

#model 2
anx.col2<- coxph(Surv(can.col.t_new, can.col.e) ~ anx.e_new + 
                   age + sex + ethn + deprivation + 
                   bmi + blm.t2d + sbp,
                 data = ds)

#Model 3
anx.col3<- coxph(Surv(can.col.t_new, can.col.e) ~ anx.e_new + 
                   age + sex + ethn + deprivation + 
                   bmi + blm.t2d + sbp +
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data = ds)

#PH assumption
phanx.col<- cox.zph(anx.col3)



### Liver Cancer ------
ds<- uk %>% 
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
                   age + sex + ethn + deprivation + 
                   bmi + blm.t2d + sbp,
                 data = ds)

#Model 3
anx.liv3<- coxph(Surv(can.liv.t_new, can.liv.e) ~ anx.e_new + 
                   age + sex + ethn + deprivation + 
                   bmi + blm.t2d + sbp +
                   smok + alc.ut + met.tot + slp.durat + ffq.pmi,
                 data = ds)

#PH assumption
phanx.liv<- cox.zph(anx.liv3)