library(broom)
library(dplyr)
library(purrr)
library(EValue)


# Overall cancer risk -----------------------------------------------------
# Put your model objects into a named list
overall_models <- list(
  Depression = depcox3,
  Anxiety = anxcox3,
  Bipolar = bplcox3,
  PTSD = ptscox3,
  Schiz = scicox3
)

# Extract HR and CI for the first coefficient of each model
overall_eval_final  <- map_df(overall_models, ~tidy(.x, exponentiate = TRUE, conf.int = TRUE) %>% 
                            slice(1) %>% 
                            select(estimate, conf.low, conf.high), 
                          .id = "MHC_Condition")  %>%
  rowwise() %>%
  mutate(
    # Create a temporary object for the E-value results
    ev_obj = list(evalues.HR(est = estimate, lo = conf.low, hi = conf.high, rare = TRUE, true = 1)),
    
    # Extract the E-value for the point estimate
    E_Value_Est = ev_obj[2, 1],
    
    # Extract the E-value for the confidence interval limit closest to the null
    E_Value_CI = pmax(ev_obj[2, 2], ev_obj[2, 3], na.rm = TRUE)) %>%
  select(-ev_obj)

print(overall_eval_final)
# Site-specific cancer risk -----------------------------------------------

site_models <- list(
  # Depression Site-Specific
  "Depression: Blood" = dep.blo3,
  "Depression: Breast" = dep.bre3,
  "Depression: Colorectal" = dep.col3,
  "Depression: Liver" = dep.liv3,
  "Depression: Lung" = dep.lun3,
  "Depression: Ovarian" = dep.ova3,
  "Depression: Prostate" = dep.pro3,
  "Depression: Uterine" = dep.ute3,
  
  # Anxiety Site-Specific
  "Anxiety: Blood" = anx.blo3,
  "Anxiety: Breast" = anx.bre3,
  "Anxiety: Colorectal" = anx.col3,
  "Anxiety: Liver" = anx.liv3,
  "Anxiety: Lung" = anx.lun3,
  "Anxiety: Ovarian" = anx.ova3,
  "Anxiety: Prostate" = anx.pro3,
  "Anxiety: Uterine" = anx.ute3
)

# Extract results and calculate E-values
site_eval_final <- map_df(site_models, ~tidy(.x, exponentiate = TRUE, conf.int = TRUE) %>% 
                               slice(1) %>% 
                               select(estimate, conf.low, conf.high), 
                             .id = "Model_Site") %>%
  rowwise() %>%
  mutate(
    # Calculate E-values (rare = TRUE because site-specific cancers are even rarer than overall)
    ev_obj = list(evalues.HR(est = estimate, lo = conf.low, hi = conf.high, rare = TRUE, true = 1)),
    E_Value_Est = ev_obj[2, 1],
    E_Value_CI = pmax(ev_obj[2, 2], ev_obj[2, 3], na.rm = TRUE)) %>%
  select(-ev_obj)


print(site_results_final)





# Table -------------------------------------------------------------------

library(gt)
library(dplyr)

# Standardize names and add a grouping column
overall_tab <- overall_eval_final %>%
  rename(Label = MHC_Condition) %>%
  mutate(Group = "Overall Cancer Risk")

site_tab <- site_eval_final %>%
  rename(Label = Model_Site) %>%
  mutate(Group = "Site-Specific Cancer Risk")

# Combine datasets
merged_results <- bind_rows(overall_tab, site_tab) %>%
  mutate(
    # Create a nice HR (95% CI) string for the table
    HR_CI = sprintf("%.2f (%.2f, %.2f)", estimate, conf.low, conf.high),
    # Round E-values to 2 decimal places
    E_Value_Est = round(E_Value_Est, 2),
    E_Value_CI = round(E_Value_CI, 2)
  )


results_table <- merged_results %>%
  select(Group, Label, HR_CI, E_Value_Est, E_Value_CI) %>%
  gt(groupname_col = "Group") %>%
  tab_header(
    title = "Table X: Sensitivity Analysis for Unmeasured Confounding showing the hazard Ratios and Corresponding E-Values for Mental Health Conditions and Cancer Risk"
  ) %>%
  cols_label(
    Label = "Condition/Site",
    HR_CI = "HR (95% CI)",
    E_Value_Est = "E-Value (Estimate)",
    E_Value_CI = "E-Value (95% CI Limit)"
  ) %>%
  # Add footnotes to explain E-values to the reader/reviewer
  tab_footnote(
    footnote = "The E-value is the minimum strength of association an unmeasured confounder would need with both the exposure and outcome to move the HR to the null.",
    locations = cells_column_labels(columns = E_Value_Est)
  ) %>%
  tab_options(
    row_group.font.weight = "bold",
    table.width = pct(100),
    column_labels.font.weight = "bold"
  )

# View the table
results_table


gtsave(results_table, filename = "output/evalue.docx")
