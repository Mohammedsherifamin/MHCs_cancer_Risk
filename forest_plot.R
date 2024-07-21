## Forest Plot of Results 
# Loading packages 
library(forestploter)

## MHCs and overall cancer risk forest plot ----

### Creating a function to extract key statistical results from the cox models
extract_cox_data <- function(model, coef_name) {
  summary_model <- summary(model)
  est <- exp(coef(model)[coef_name])
  ci <- exp(confint(model)[coef_name, ])
  p_value <- summary_model$coefficients[coef_name, "Pr(>|z|)"]
  nevent <- model$nevent
  n <- model$n
  data.frame(
    Model = coef_name,
    est = est,
    low = ci[1],
    hi = ci[2],
    p_value = p_value,
    nevent = nevent,
    n = n
  )
}

# Specifying the coefficients to be extracted
models <- list(depcox1 = "dep.e_new1", depcox2 = "dep.e_new1", depcox3 = "dep.e_new1",
               anxcox1 = "anx.e_new1", anxcox2 = "anx.e_new1", anxcox3 = "anx.e_new1",
               bplcox1 = "bpl.e_new1", bplcox2 = "bpl.e_new1", bplcox3 = "bpl.e_new1",
               scicox1 = "sci.e_new1", scicox2 = "sci.e_new1", scicox3 = "sci.e_new1",
               ptscox1 = "srd.pts.e_new1", ptscox2 = "srd.pts.e_new1", ptscox3 = "srd.pts.e_new1")

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
                          ifelse(p_value < 0.05, paste0(round(p_value, 4), "*"), 
                                 round(p_value, 4))))

model_names <- c("Model 1", "Model 2", "Model 3")

# Calculate the indices for each row
indices <- seq_len(nrow(results)) - 1 

# Use modulo to cycle through model names
results$Model <- model_names[(indices %% 3) + 1]

results <- results %>%
  mutate(Model = paste0("       ", Model))

# Create the new rows with the specified names and empty cells
new_rows <- data.frame(Model = c("Depressive Disorders", "", "", "", "Anxiety disorders", "", "", "", "Bipolar Disorders", "", "", "", "Schizophrenia", "", "", "", "PTSD"),
                       matrix("", nrow = 17, ncol = ncol(results) - 1))

# Set the column names to match those of the existing data frame
colnames(new_rows) <- colnames(results)

# Insert the new rows 
results <- rbind(new_rows[1,], results)
results <- rbind(results[1:4,], new_rows[5,], results[5:nrow(results),])
results <- rbind(results[1:8,], new_rows[9,], results[9:nrow(results),])
results <- rbind(results[1:12,], new_rows[13,], results[13:nrow(results),])
results <- rbind(results[1:16,], new_rows[17,], results[17:nrow(results),])

results$est<- as.numeric(results$est)
results$low<- as.numeric(results$low)
results$hi<- as.numeric(results$hi)

# Adjust headings name
names(results)[names(results) == "nevent"] <- "Number of Events"
names(results)[names(results) == "n"] <- "Effective Sample Size"

# Plot the dataframe
f_mhc_all <- forest(
  results[, c("Model",
              "Number of Events" ," ", "HR (95% CI)", "p_value")],
  est = results$est,
  lower = results$low, 
  upper = results$hi,
  ci_column = 3,
  ref_line = 1,
  arrow_lab = c("Lower Cancer Risk", "Higher Cancer Risk"),
  xlim = c(0.8, 2.5),
  ticks_at = c(0.8,1, 1.5, 2, 2.5),
  xlab = "HR",
  title = "Forest Plot of Cox Models for MHCs and Overall Cancer"
)

### Adusting the plot
# Adding borders
f_mhc_all <-add_border(f_mhc_all, part = "header", row = 1, 
                       gp =gpar(lwd= 1))

#Adding footnotes 
f_mhc_all <- add_text(
  f_mhc_all, 
  text = "
HR: Hazard Ratio; CI: Confidence Interval; PTSD: Post-Traumatic Stress Disorders.\nModel 1 adjusted for sociodemographic factors; Model 2 adjusted for sociodemographic  and lifestyle factors;\nModel 3 adjusted for sociodemographic factors, lifestyle and health-related factors.\n*Significance level set at alpha = 0.05.",
  row = nrow(results) + 3, 
  col = 1:3, 
  just = "left",
  gp = gpar(fontsize = 7.8, col = "black", fontface= 3),
  padding = unit(4, "mm")
)

### Save object for dissertation
ggplot2::ggsave(filename = "pic/f_mhc_all.png", plot = f_mhc_all,
                dpi = 300,
                width = 6.8, height = 6.8, units = "in")

## Sex Stratified Forest Plot -----
#Creating function to extract data from the sex stratified cox models
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

# Specifying the coefficents to be extracted
models <- list(fdepcox = "dep.e_new1", mdepcox = "dep.e_new1",
               fbplcox = "bpl.e_new1", mbplcox = "bpl.e_new1")

results_list <- lapply(names(models), function(model_name) {
  model <- get(model_name)
  coef_name <- models[[model_name]]
  extract_cox_data(model, coef_name)
})

results <- do.call(rbind, results_list)

results$` ` <- paste(rep(" ", 20), collapse = " ")

# Create a CI to display
results$`HR (95% CI)` <- sprintf("%.2f (%.2f - %.2f)", results$est, 
                                 results$low, results$hi)

# adjust p-value format 
results <- results %>%
  mutate(p_value = ifelse(p_value < 0.001, "<0.001*", 
                          ifelse(p_value < 0.05, paste0(round(p_value, 4), "*"), 
                                 round(p_value, 4))))

#Adjust models names 
model_names <- c("Female ", "Male")

# Calculate the indices for each row
indices <- seq_len(nrow(results)) - 1  

# Use modulo to cycle through model names
results$Model <- model_names[(indices %% 2) + 1]
results <- results %>%
  mutate(Model = paste0("       ", Model))

# Create the new rows with the specified names and empty cells to act as headings 
new_rows <- data.frame(Model = c("Depressive Disorders", "", "", "Bipolar Disorders"),
                       matrix("", nrow = 4, ncol = ncol(results) - 1))

# Set the column names to match those of the existing data frame
colnames(new_rows) <- colnames(results)

# Insert the new rows at the specified positions
results <- rbind(new_rows[1,], results)
results <- rbind(results[1:3,], new_rows[4,], results[4:nrow(results),])
results <- rbind(results[1:6,], new_rows[7,], results[7:nrow(results),])
results <- rbind(results[1:9,], new_rows[10,], results[10:nrow(results),])
results <- rbind(results[1:12,], new_rows[13,], results[13:nrow(results),])

results$est<- as.numeric(results$est)
results$low<- as.numeric(results$low)
results$hi<- as.numeric(results$hi)

# Adjust headings names 
names(results)[names(results) == "nevent"] <- "Number of Events"
names(results)[names(results) == "n"] <- "Total Observations"

# plot the data frame 
f_mhc_sex <- forest(
  results[, c("Model","Number of Events", " ", "HR (95% CI)", "p_value")],
  est = results$est,
  lower = results$low, 
  upper = results$hi,
  ci_column = 3,
  ref_line = 1,
  arrow_lab = c("Lower Cancer Risk", "Higher Cancer Risk"),
  xlim = c(0.5, 2),
  ticks_at = c(0.5, 1, 1.5, 2),
  xlab = "HR",
  title = "Forest Plot for MHCs and Overall Cancer Stratified by Sex",
)

# Adjust the plot 
## Adding borders
f_mhc_sex <-add_border(f_mhc_sex, part = "header", row = 1, 
                       gp =gpar(lwd= 1))

## Adding footnotes
f_mhc_sex  <- add_text(
  f_mhc_sex, 
  text = "
HR: Hazard Ratio; CI: Confidence Interval\nAll Model are adjusted for sociodemographic factors, lifestyle and health-related factors.\n*Significance level set at alpha = 0.05",  row = nrow(results) + 3, # Adjust row number based on your plot
  col = 1:3, 
  just = "left",
  gp = gpar(fontsize = 9, col = "black", fontface= 3), 
  padding = unit(4, "mm")
)

## Saving object for dissertation
ggplot2::ggsave(filename = "pic/mhc_forest_sex.png", plot = f_mhc_sex,
                dpi = 300,
                width = 7.5, height = 4, units = "in")



# Depression and Site specific cancer -----
#Creating function to extract data from depression cox models
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

# Specifying the coefficients to be extracted
models <- list(dep.bre1 = "dep.e_new1", dep.bre2 = "dep.e_new1", dep.bre3= "dep.e_new1",
               dep.ova1 = "dep.e_new1", dep.ova2 = "dep.e_new1", dep.ova3= "dep.e_new1",
               dep.ute1 = "dep.e_new1", dep.ute2 = "dep.e_new1", dep.ute3= "dep.e_new1",
               dep.pro1 = "dep.e_new1", dep.pro2 = "dep.e_new1", dep.pro3= "dep.e_new1",
               dep.lun1 = "dep.e_new1", dep.lun2 = "dep.e_new1", dep.lun3= "dep.e_new1",
               dep.blo1 = "dep.e_new1", dep.blo2 = "dep.e_new1", dep.blo3= "dep.e_new1",
               dep.col1 = "dep.e_new1", dep.col2 = "dep.e_new1", dep.col3= "dep.e_new1",
               dep.liv1 = "dep.e_new1", dep.liv2 = "dep.e_new1", dep.liv3= "dep.e_new1")

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
results$`HR (95% CI)` <- sprintf("%.2f (%.2f - %.2f)", results$est, results$low,
                                 results$hi)

# Change p-value format
results <- results %>%
  mutate(p_value = ifelse(p_value < 0.001, "<0.001*", 
                          ifelse(p_value < 0.05, paste0(round(p_value, 4), "*"), 
                                 round(p_value, 4))))

model_names <- c("Model 1", "Model 2", "Model 3")

# Calculate the indices for each row
indices <- seq_len(nrow(results)) - 1 

# Use modulo to cycle through model names
results$Model <- model_names[(indices %% 3) + 1]
results <- results %>%
  mutate(Model = paste0("       ", Model))

# Create the new rows with the specified names and empty cells
new_rows <- data.frame(Model = c("Breast Cancer", "", "", "", "Ovarian Cancer", "", "", "", "Uterine Cancer", "", "", "", "Prostate Cancer", "", "", "", "Lung Cancer",
                                 "", "", "", "Blood Cancer", "", "", "", "Colorectal Cancer", "", "", "","Liver Cancer"),
                       matrix("", nrow = 29, ncol = ncol(results) - 1))

# Set the column names to match those of the existing data frame
colnames(new_rows) <- colnames(results)

# Insert the new rows at the specified positions
results <- rbind(new_rows[1,], results)
results <- rbind(results[1:4,], new_rows[5,], results[5:nrow(results),])
results <- rbind(results[1:8,], new_rows[9,], results[9:nrow(results),])
results <- rbind(results[1:12,], new_rows[13,], results[13:nrow(results),])
results <- rbind(results[1:16,], new_rows[17,], results[17:nrow(results),])
results <- rbind(results[1:20,], new_rows[21,], results[21:nrow(results),])
results <- rbind(results[1:24,], new_rows[25,], results[25:nrow(results),])
results <- rbind(results[1:28,], new_rows[29,], results[29:nrow(results),])

results$est<- as.numeric(results$est)
results$low<- as.numeric(results$low)
results$hi<- as.numeric(results$hi)

names(results)[names(results) == "nevent"] <- "Number of Events"

# plot the dataframe
f_dep <- forest(
  results[, c("Model", "Number of Events", " ", "HR (95% CI)", "p_value")],
  est = results$est,
  lower = results$low, 
  upper = results$hi,
  ci_column = 3,
  ref_line = 1,
  arrow_lab = c("Lower Cancer Risk", "Higher Cancer Risk"),
  xlim = c(0, 3.5),
  ticks_at = c(0.5, 1, 1.5, 2,3, 3.5),
  xlab = "HR",
  title = "Forest Plot of Cox Models for Depressive Diosrders and site-specific Cancer",
)

## Adjusting the plot 
# Adding borders 
f_dep <-add_border(f_dep,  part = "header", row = 1, 
                   gp =gpar(lwd= 1))

# Adding footnotes
f_dep<- add_text(
  f_dep, 
  text = "
HR: Hazard Ratio; CI: Confidence Interval.\nModel 1 adjusted for sociodemographic factors; Model 2 adjusted for sociodemographic  and lifestyle factors;\nModel 3 adjusted for sociodemographic factors, lifestyle and health-related factors.\n*Significance level set at alpha = 0.05.",
  row = nrow(results) + 3, 
  col = 1:3,
  just = "left",
  gp = gpar(fontsize = 9, col = "black", fontface= 3), 
  padding = unit(4, "mm")
)

# Saving object
ggplot2::ggsave(filename = "pic/f_dep.png", plot = f_dep,
                dpi = 300,
                width = 8.5, height = 10, units = "in")

# Anxiety and Site Specific cancer --------
#Creating function to extract data from anxiety cox models
xtract_cox_data <- function(model, coef_name) {
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

# Specifying coefficient to be plotted 
models <- list(anx.bre1 = "anx.e_new1", anx.bre2 = "anx.e_new1", anx.bre3= "anx.e_new1",
               anx.ova1 = "anx.e_new1", anx.ova2 = "anx.e_new1", anx.ova3= "anx.e_new1",
               anx.ute1 = "anx.e_new1", anx.ute2 = "anx.e_new1", anx.ute3= "anx.e_new1",
               anx.pro1 = "anx.e_new1", anx.pro2 = "anx.e_new1", anx.pro3= "anx.e_new1",
               anx.lun1 = "anx.e_new1", anx.lun2 = "anx.e_new1", anx.lun3= "anx.e_new1",
               anx.blo1 = "anx.e_new1", anx.blo2 = "anx.e_new1", anx.blo3= "anx.e_new1",
               anx.col1 = "anx.e_new1", anx.col2 = "anx.e_new1", anx.col3= "anx.e_new1",
               anx.liv1 = "anx.e_new1", anx.liv2 = "anx.e_new1", anx.liv3= "anx.e_new1")

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
                          ifelse(p_value < 0.05, paste0(round(p_value, 4), "*"), 
                                 round(p_value, 4))))

model_names <- c("Model 1", "Model 2", "Model 3")

# Calculate the indices for each row
indices <- seq_len(nrow(results)) - 1 

# Use modulo to cycle through model names
results$Model <- model_names[(indices %% 3) + 1]
results <- results %>%
  mutate(Model = paste0("       ", Model))

# Create the new rows with the specified names and empty cells
new_rows <- data.frame(Model = c("Breast Cancer", "", "", "", "Ovarian Cancer", "", "", "", "Uterine Cancer", "", "", "", "Prostate Cancer", "", "", "", "Lung Cancer",
                                 "", "", "", "Blood Cancer", "", "", "", "Colorectal Cancer", "", "", "","Liver Cancer"),
                       matrix("", nrow = 29, ncol = ncol(results) - 1))

# Set the column names to match those of the existing data frame
colnames(new_rows) <- colnames(results)

# Insert the new rows at the specified positions
results <- rbind(new_rows[1,], results)
results <- rbind(results[1:4,], new_rows[5,], results[5:nrow(results),])
results <- rbind(results[1:8,], new_rows[9,], results[9:nrow(results),])
results <- rbind(results[1:12,], new_rows[13,], results[13:nrow(results),])
results <- rbind(results[1:16,], new_rows[17,], results[17:nrow(results),])
results <- rbind(results[1:20,], new_rows[21,], results[21:nrow(results),])
results <- rbind(results[1:24,], new_rows[25,], results[25:nrow(results),])
results <- rbind(results[1:28,], new_rows[29,], results[29:nrow(results),])

results$est<- as.numeric(results$est)
results$low<- as.numeric(results$low)
results$hi<- as.numeric(results$hi)

# Adusting heading names
names(results)[names(results) == "nevent"] <- "Number of Events"

# plotting the dataframe
f_anx <- forest(
  results[, c("Model", "Number of Events", " ", "HR (95% CI)", "p_value")],
  est = results$est,
  lower = results$low, 
  upper = results$hi,
  ci_column = 3,
  ref_line = 1,
  arrow_lab = c("Lower Cancer Risk", "Higher Cancer Risk"),
  xlim = c(0, 4),
  ticks_at = c(0.5, 1, 2,3, 4),
  xlab = "HR",
  title = "Forest Plot of Cox Models for Anxiety Diosrders and site-specific Cancer",
)

# Adjusting the plot 
## Adding borders
f_anx <-add_border(f_anx,  part = "header", row = 1, 
                   gp =gpar(lwd= 1))

# adding footnotes
f_anx<- add_text(
  f_anx, 
  text = "
HR: Hazard Ratio; CI: Confidence Interval.\nModel 1 adjusted for sociodemographic factors; Model 2 adjusted for sociodemographic  and lifestyle factors;\nModel 3 adjusted for sociodemographic factors, lifestyle and health-related factors.\n*Significance level set at alpha = 0.05.",
  row = nrow(results) + 3, 
  col = 1:3, 
  just = "left", 
  gp = gpar(fontsize = 9, col = "black", fontface= 3), 
  padding = unit(4, "mm")
)

# saving object 
ggplot2::ggsave(filename = "pic/f_anx.png", plot = f_anx,
                dpi = 250,
                width = 8, height = 10, units = "in")



# Bipolar Disorders and site specific cancer -------------
# Creating function to extract data from Bipolar disorders cox models 
xtract_cox_data <- function(model, coef_name) {
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

# Specifying the coefficients to be plotted 
models <- list(bpl.bre1 = "bpl.e_new1", bpl.bre2 = "bpl.e_new1", bpl.bre3= "bpl.e_new1",
               bpl.pro1 = "bpl.e_new1", bpl.pro2 = "bpl.e_new1", bpl.pro3= "bpl.e_new1",
               bpl.blo1 = "bpl.e_new1", bpl.blo2 = "bpl.e_new1", bpl.blo3= "bpl.e_new1")

# Extract data from all models
results_list <- lapply(names(models), function(model_name) {
  model <- get(model_name)
  coef_name <- models[[model_name]]
  extract_cox_data(model, coef_name)
})

# Combine all results into a data frame
results <- do.call(rbind, results_list)

# Add a blank column for the forest plot to display CI
results$` ` <- paste(rep(" ", 20), collapse = " ")

# Create a confidence interval column to display
results$`HR (95% CI)` <- sprintf("%.2f (%.2f - %.2f)", results$est, results$low, results$hi)

# Change p-value format
results <- results %>%
  mutate(p_value = ifelse(p_value < 0.001, "<0.001*", 
                          ifelse(p_value < 0.05, paste0(round(p_value, 4), "*"), 
                                 round(p_value, 4))))

# Adding model names
model_names <- c("Model 1", "Model 2", "Model 3")

# Calculate the indices for each row
indices <- seq_len(nrow(results)) - 1  

# Use modulo to cycle through model names
results$Model <- model_names[(indices %% 3) + 1]
results <- results %>%
  mutate(Model = paste0("       ", Model))

# Create the new rows with the specified names and empty cells
new_rows <- data.frame(Model = c("Breast Cancer", "", "", "", "Prostate Cancer", "", "", "", "Blood Cancer"),
                       matrix("", nrow = 9, ncol = ncol(results) - 1))

# Set the column names to match those of the existing data frame
colnames(new_rows) <- colnames(results)

# Insert the new rows at the specified positions
results <- rbind(new_rows[1,], results)
results <- rbind(results[1:4,], new_rows[5,], results[5:nrow(results),])
results <- rbind(results[1:8,], new_rows[9,], results[9:nrow(results),])
results <- rbind(results[1:12,], new_rows[13,], results[13:nrow(results),])
results <- rbind(results[1:16,], new_rows[17,], results[17:nrow(results),])
results <- rbind(results[1:20,], new_rows[21,], results[21:nrow(results),])
results <- rbind(results[1:24,], new_rows[25,], results[25:nrow(results),])
results <- rbind(results[1:28,], new_rows[29,], results[29:nrow(results),])

results$est<- as.numeric(results$est)
results$low<- as.numeric(results$low)
results$hi<- as.numeric(results$hi)

# Adjusting heading name 
names(results)[names(results) == "nevent"] <- "Number of Events"

#Plotting the data frame
f_bpl <- forest(
  results[, c("Model", "Number of Events", " ", "HR (95% CI)", "p_value")],
  est = results$est,
  lower = results$low, 
  upper = results$hi,
  ci_column = 3,
  ref_line = 1,
  arrow_lab = c("Lower Cancer Risk", "Higher Cancer Risk"),
  xlim = c(0, 3),
  ticks_at = c(0.5, 1, 2,3),
  xlab = "HR",
  title = "Forest Plot of Cox Models for Bipolar Diosrders and site-specific Cancer",
)

# adjusting the plot 
## adding borders 
f_bpl <-add_border(f_bpl,  part = "header", row = 1, 
                   gp =gpar(lwd= 1))

# Adding footnotes 
f_bpl<- add_text(
  f_bpl, 
  text = "
HR: Hazard Ratio; CI: Confidence Interval.\nModel 1 adjusted for sociodemographic factors; Model 2 adjusted for sociodemographic  and lifestyle factors;\nModel 3 adjusted for sociodemographic factors, lifestyle and health-related factors.\n*Significance level set at alpha = 0.05.",
  row = nrow(results) + 3,
  col = 1:3, 
  just = "left",
  gp = gpar(fontsize = 9, col = "black", fontface= 3),
  padding = unit(4, "mm")
)

# Saving the object
ggplot2::ggsave(filename = "pic/f_bpl.png", plot = f_bpl,
                dpi = 250,
                width = 8, height = 5.5, units = "in")
