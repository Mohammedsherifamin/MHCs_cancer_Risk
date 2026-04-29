rm(list = ls())
source("scripts/1.data_prepration.R")

#Loading Packages 
library(gtsummary)
library(flextable)
library(tidyverse)
library(gt)
library(broom)
library(naniar)

## Check the amount of missing data 
missing_counts <- colSums(is.na(uk))   

print(missing_counts)

#select only relevant variables for the main analysis
uk <- uk %>%
  dplyr::select(eid, age, sex,ethn,
                deprivation, smok,alc.ut, slp.durat, 
                met.tot, ffq.pmi, bmi, blm.t2d, sbp, menopause, mhc.e,
                dep.e, anx.e, bpl.e, sci.e, srd.pts.e, dep.t, anx.t, bpl.t, 
                sci.t, srd.pts.t, can.e, can.t, can.bre.e, can.ova.e, can.ute.e, 
                can.pro.e, can.lun.e, can.blo.e, can.col.e, can.liv.e, can.bre.t, 
                can.ova.t, can.ute.t, can.pro.t, can.lun.t ,can.blo.t , 
                can.col.t, can.liv.t)

nrow(uk)

# --------------------------- Missing data -------------------------------
# create indicator for included vs excluded due to missing data
uk_compare <- uk %>%
  mutate(
    inclusion_status = if_else(
      complete.cases(across(all_of(everything(.)))),
      "Included in analysis",
      "Excluded due to missing data"
    ),
    mhc.e = case_when(
      mhc.e == 0 ~ "Undiagnosed",
      mhc.e == 1 ~ "Diagnosed",
      TRUE ~ NA_character_
    )
  )

# check numbers
table(uk_compare$inclusion_status)


# Creating a table for included vs excluded
base_char <-  c(
  "age",
  "sex",
  "ethn",
  "deprivation",
  "smok",
  "alc.ut" ,
  "slp.durat",
  "met.tot",
  "ffq.pmi",
  "bmi" ,
  "blm.t2d",
  "sbp"    ,
  "menopause",
  "can.t", 
  "mhc.e"
)


tbl_included_excluded <- uk_compare %>%
  tbl_summary(
    by = inclusion_status,
    include = all_of(base_char),
    missing = "no",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      c("alc.ut", "can.t") ~ "{median} ({p25}, {p75})"
    ),
    label = list(
      age ~ "Age (years)",
      sex ~ "Sex",
      ethn ~ "Ethnicity",
      deprivation ~ "Deprivation",
      smok ~ "Smoking Status",
      alc.ut ~ "Alcohol (units/week)",
      slp.durat ~ "Sleep Duration (hour)",
      met.tot ~ "Total physical activity (MET-min/week)",
      ffq.pmi ~ "Processed Meat Consumption",
      bmi ~ "BMI (kg/m2)",
      blm.t2d ~ "Self-reported Diabetes",
      sbp ~ "SBP (mmHg)",
      menopause ~ "Menopausal Status",
      mhc.e ~ "Mental health condition diagnosis",
      can.t ~ "Follow-up Time (year)"),
    
    type = list(
      age ~ "continuous",
      deprivation ~ "continuous",
      slp.durat ~ "continuous",
      met.tot ~ "continuous",
      bmi ~ "continuous",
      sbp ~ "continuous",
      alc.ut ~ "continuous",
      can.t ~ "continuous")) %>%
  
  add_overall() %>%
  
  add_p(
    test = list(
      age ~ "t.test",
      deprivation ~ "t.test",
      slp.durat ~ "t.test",
      met.tot ~ "t.test",
      bmi ~ "t.test",
      sbp ~ "t.test",
      alc.ut ~ "wilcox.test",
      can.t ~ "wilcox.test",
      all_categorical() ~ "chisq.test")) 

# export to Word
tbl_included_excluded %>%
  as_flex_table() %>%
  save_as_docx(path = "output/included_vs_excluded_table.docx")

tbl_included_excluded


naniar::miss_var_pct(uk)
miss_var_summary(uk) %>% gt()
## ---------------------------- Final analytic sample -----------------------------
# Final analytic sampel baseline characterestics table

## Due to the small number of missing data, we decided to drop missing data and preform a 
#complete cases analysis
uk <- uk %>%  drop_na()

nrow(uk)

## Table summary of baseline characteristics -----
base_char <-  c(
  "age",
  "sex",
  "ethn",
  "deprivation",
  "smok",
  "alc.ut" ,
  "slp.durat",
  "met.tot",
  "ffq.pmi",
  "bmi" ,
  "blm.t2d",
  "sbp"    ,
  "menopause",
  "can.t", 
  "mhc.e"
)

tbl_sum <- tbl_summary(
  uk %>% mutate(mhc.e = case_when(
    mhc.e==0 ~ "Undiagnosed", 
    mhc.e==1 ~ "Diagnosed")),
  by = mhc.e,
  missing = "no",
  include = setdiff(base_char, "mhc.e"),
  statistic = list(
    all_continuous() ~ "{mean} ({sd})",
    c("can.t", "alc.ut") ~ "{median} ({p25}, {p75})"),
  type = list(c(slp.durat) ~ "continuous"),
  label = list(
    age ~ "Age (years)",
    sex ~ "Sex",
    ethn ~ "Ethnicity",
    deprivation ~ "Deprivation",
    smok ~ "Smoking Status",
    alc.ut ~ "Alcohol (units/week)",
    slp.durat ~ "Sleep Duration (hour)",
    met.tot ~ "Total physical activity (MET-min/week)",
    ffq.pmi ~ "Processed Meat Consumption",
    bmi ~ "BMI (kg/m)",
    blm.t2d ~ "Self-reported Diabetes",
    sbp ~ "SBP (mmHg)",
    menopause ~ "Menopausal Status",
    can.t ~ "Follow-up Time (year)")) %>%
  add_overall() %>%
  add_p(
    test = list(
      age ~ "t.test",
      deprivation ~ "t.test",
      slp.durat ~ "t.test",
      met.tot ~ "t.test",
      bmi ~ "t.test",
      sbp ~ "t.test",
      alc.ut ~ "wilcox.test",
      can.t ~ "wilcox.test",
      all_categorical() ~ "chisq.test"))

tbl_sum <- tbl_sum %>% as_flex_table() 
tbl_sum %>% save_as_docx(path = "output/summary_overall.docx")

#Stratifying baseline characteristics by sex------
tbl_sum.sex <- tbl_summary(
  uk,
  by = sex,
  missing = "no",
  include = setdiff(base_char, "mhc.e"),
  statistic = list(
    all_continuous() ~ "{mean} ({sd})",
    c("can.t", "alc.ut") ~ "{median} ({p25}, {p75})"),
  type = list(c(slp.durat) ~ "continuous"),
  label = list(
    age ~ "Age (years)",
    ethn ~ "Ethnicity",
    deprivation ~ "Deprivation",
    smok ~ "Smoking Status",
    alc.ut ~ "Alcohol (units/week)",
    slp.durat ~ "Sleep Duration (hour)",
    met.tot ~ "Total physical activity (MET-min/week)",
    ffq.pmi ~ "Processed Meat Consumption",
    bmi ~ "BMI (kg/m)",
    blm.t2d ~ "Self-reported Diabetes",
    sbp ~ "SBP (mmHg)",
    menopause ~ "Menopausal Status",
    can.t ~ "Follow-up Time (year)")) %>%
  add_overall() %>%
  add_p(
    test = list(
      age ~ "t.test",
      deprivation ~ "t.test",
      slp.durat ~ "t.test",
      met.tot ~ "t.test",
      bmi ~ "t.test",
      sbp ~ "t.test",
      alc.ut ~ "wilcox.test",
      can.t ~ "wilcox.test",
      all_categorical() ~ "chisq.test"))

print(tbl_sum.sex)

tbl_sum.sex <- tbl_sum.sex %>% as_flex_table() 

tbl_sum.sex %>% save_as_docx(path = "output/summary_sex.docx")


