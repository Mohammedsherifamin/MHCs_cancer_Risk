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


#select only relevant variables, including phq this time.
uk <- uk %>%
  dplyr::select(eid, age, sex,ethn, phq,
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
ds <- uk %>%  drop_na()

nrow(ds)



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
                , data= subset(ds, sex=="F"))

#uterine 
sensute<- coxph(Surv(can.ute.t, can.ute.e) ~ 
                  phqs +age + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp + menopause 
                , data= subset(ds, sex=="F"))

#Prostate cancer 
senspro<- coxph(Surv(can.pro.t, can.pro.e) ~ 
                  phqs +age + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp 
                , data= subset(ds, sex=="M"))

# lung cancer 
senslun<- coxph(Surv(can.lun.t, can.lun.e) ~ phqs+ 
                  age +sex + ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp , data=ds)

#Blood Cancer 
sensblo<- coxph(Surv(can.blo.t, can.blo.e) ~ 
                  phqs +age + sex+ ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp , data= ds)

#Colon Cancer 
senscol<- coxph(Surv(can.col.t, can.col.e) ~ 
                  phqs +age + sex+  ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp , data= ds)

#Liver Cancer 
sensliv<- coxph(Surv(can.liv.t, can.liv.e) ~ 
                  phqs +age + sex+  ethn + deprivation + 
                  smok + alc.ut + met.tot + slp.durat + ffq.pmi +
                  bmi + blm.t2d + sbp , data= ds)

# sensitivity analysis forest plot ----
## Creating a function to extract data from the sensitivity analysis cox models 
extract_cox_data <- function(model, coef_name) {
  summary_model <- summary(model)
  est <- exp(coef(model)[coef_name])
  ci <- exp(confint(model)[coef_name, ])
  p_value <- summary_model$coefficients[coef_name, "Pr(>|z|)"]
  nevent <- model$nevent
  data.frame(
    Model = coef_name,
    est = est,
    low = ci[1],
    hi = ci[2],
    p_value = p_value,
    nevent = nevent
  )
}

# Specifying the coefficeints to be extracted 
models <- list(sensall = "phqs",sensbre = "phqs",sensova = "phqs",
               sensute = "phqs",senspro = "phqs",senslun = "phqs",
               sensblo = "phqs",senscol = "phqs",sensliv = "phqs")

# Extract data from all models
results_list <- lapply(names(models), function(model_name) {
  model <- get(model_name)
  coef_name <- models[[model_name]]
  extract_cox_data(model, coef_name)
})

# Combine all results into a single data frame
results <- do.call(rbind, results_list)

# Add a blank column for the forest plot to display CI
results$` ` <- paste(rep(" ", 20), collapse = " ")

# Create a confidence interval column to display
results$`HR (95% CI)` <- sprintf("%.2f (%.2f - %.2f)", results$est, results$low, results$hi)

# Change p-value format
results <- results %>%
  mutate(p_value = ifelse(p_value < 0.001, "<0.001*", 
                          ifelse(p_value < 0.01, paste0(round(p_value, 3), "*"), 
                                 ifelse(p_value < 0.05, paste0(round(p_value, 2), "*"), 
                                        round(p_value, 2)))))
# Creating new names to be displayed in the data frame
model_names <- c("Overall Cancer", "Breast Cancer", "Ovarian Cancer",  "Uterine Cancer", "Prostate Cancer", 
                 "Lung Cancer", "Blood Cancer", "Colorectal Cancer","Liver Cancer")

#Adding names to the dataframe 
results$Model <- model_names

results$est<- as.numeric(results$est)
results$low<- as.numeric(results$low)
results$hi<- as.numeric(results$hi)

# Adjsitng heading names 
names(results)[names(results) == "nevent"] <- "Number of Events"

# Ploting the dataframe 
f_sens <- forest(
  results[, c("Model", "Number of Events", " ","HR (95% CI)", "p_value")],
  est = results$est,
  lower = results$low, 
  upper = results$hi,
  ci_column = 3,
  ref_line = 1,
  arrow_lab = c("Lower Cancer Risk", "Higher Cancer Risk"),
  xlim = c(0.9, 1.2),
  ticks_at = c(0.9,1,1.1,1.2),
  xlab = "HR",
  title = "Sensitivity Analysis", 
  theme = forest_theme(base_family = "Helvetica")
)

# Adjusting the plot 
## adding borders
f_sens <-add_border(f_sens, part = "header", row = 1, 
                    gp =gpar(lwd= 1, fontfamily = "Helvetica"))

# Adding footnotes
f_sens  <- add_text(
  f_sens , 
  text = "
HR: Hazard Ratio; CI: Confidence Interval; \nThe standardised 4-item Patient Health Questionnaire (PHQ-4) scores represented the main independent \nvariable across all models; All Model are adjusted for sociodemographic, health-related, and lifestyle factors.\n*Significance level set at alpha = 0.05",  row = nrow(results) + 3, # Adjust row number based on your plot
  col = 1:3, 
  just = "left",
  gp = gpar(fontsize = 9, col = "black", fontface= 3, fontfamily = "Helvetica"),
  padding = unit(4, "mm")
)

# Saving the object 
ggplot2::ggsave(filename = "output/new/Supplementary Figure 2-new.pdf", plot = f_sens,
                width = 8, height = 5, units = "in")
