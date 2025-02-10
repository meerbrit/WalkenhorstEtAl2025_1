### Calculate milestone ages for FLUT study subjects (by ID!) and include all predictors ####
### BWalkenhorst, 2024 #######################################################################
#### SETUP ####
# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #includes hidden objects.
gc() 

# load all necessary libraries
library(ggplot2) 
library(ggpubr) 
library(tidyverse)
library(tidybayes) 
library(brms)
library(bayestestR) #e.g. diagnostic_posterior
library(bayesplot)
library(ggokabeito) # colour palette
library(emmeans) # 
library(extrafont)
# font_import() # use on first use

set.seed(23)

PROP_data <- read_rds('../data/PROP_data.rds')
B_prop <- readRDS('../models/B_prop.rds')

MIN_AGE = 1
MAX_AGE = 180
PROP_DIFF = 0.01
NUM_SAMPLES_PER_ID = 100

# Custom ggplot theme 
theme_clean <- function() {
  theme_minimal(base_family='Calibri') +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = rel(2), hjust = 0.5),
      axis.title = element_text(face = "bold", size = rel(1.75)),
      axis.text = element_text(face = "bold", size= rel(1.25)),
      strip.text = element_text(face = "bold", size = rel(1.5), color='white'),
      strip.background = element_rect(fill = "grey80", color = NA),
      legend.title = element_text(face = "bold", size = rel(1.25)),
      legend.text = element_text(face = 'italic', size = rel(1)))
}

# Functions for extracting posterior proportions and determine milestone ages
extract_posterior_call_props <- function(row) {
  rec_age_c <- seq(MIN_AGE, MAX_AGE, by = 1)
  age_z_vals <- (rec_age_c - mean(PROP_data$REC_AGE)) / sd(PROP_data$REC_AGE)
  
  # Expand AGE_z while keeping all other predictors constant
  new_data <- expand_grid(
    AGE_z = age_z_vals
  ) %>%
    mutate(
      REC_AGE = rec_age_c,  # Assign corresponding REC_AGE
      TREATMENT = row$TREATMENT,
      SEX = row$SEX,
      WEIGHT_z = row$WEIGHT_z,
      COMP_NORM_z = row$COMP_z,
      GS_z = row$GS_z,
      RAIN_z = row$RAIN_z,
      MOTHER_ID = row$MOTHER_ID,
      LITTER_CODE = row$LITTER_CODE,
      ID = row$ID,
      Total_calls = 1,
      .row = row_number()  # Assign `.row` index for merging
    )
  
  # Generate posterior predictions
  PROP_pred <- B_prop %>%
    epred_draws(newdata = new_data, re_formula = NULL, allow_new_levels = FALSE, robust = TRUE, ndraws = NULL)  # Keep ndraws draws
  
  # Merge posterior predictions with data
  PROP_pred_tidy <- PROP_pred %>%
    left_join(new_data, by = ".row") %>%  
    mutate(
      Call_type = .category,  # Preserve call type
      Call_prop = .epred  # Rename posterior prediction column
    ) %>%
    select(-.epred, -.category, -.row)  # Remove redundant columns
  
  # Remove duplicated columns
  PROP_pred_tidy <- PROP_pred_tidy %>%
    select(!ends_with(".y")) %>%
    rename_with(~ str_remove(., "\\.x$"))
  
  return(PROP_pred_tidy)
}

find_transition_ages <- function(call_data, ID_DATA, num_samples = NUM_SAMPLES_PER_ID) {
  call_data <- call_data %>%
    mutate(Call_type = case_when(
      .category %in% c("SumDIG") ~ "DIG",
      .category %in% c("SumBEG") ~ "REP",
      .category %in% c("SumCC") ~ "CC",
      TRUE ~ as.character(.category)  
    )) %>%
    select(-.category)  # Remove redundant column after renaming
  
  # Group & reshape to wide format
  call_data <- call_data %>%
    group_by(ID, REC_AGE, Call_type, .draw) %>%
    summarise(Call_prop = mean(Call_prop), .groups = "drop") %>% 
    pivot_wider(names_from = Call_type, values_from = Call_prop) %>%
    na.omit()  
  
  # Normalize proportions per age so they sum to 1
  call_data <- call_data %>%
    mutate(TOTAL_PROP = DIG + REP + CC) %>%
    mutate(
      DIG = DIG / TOTAL_PROP,
      REP = REP / TOTAL_PROP,
      CC = CC / TOTAL_PROP
    ) %>%
    select(-TOTAL_PROP)
  
  # Ensure DIG > REP with a minimum threshold and enforce age constraints
  semi_age <- call_data %>%
    filter(DIG > REP + PROP_DIFF) %>%
    group_by(.draw) %>%
    summarise(SEMI_AGE = min(REC_AGE, na.rm = TRUE), .groups = "drop") %>%
    mutate(SEMI_AGE = ifelse(SEMI_AGE < 30, NA, SEMI_AGE))
  
  # Find the peak DIG age (max DIG proportion) after SEMI_AGE
  peak_dig_age <- call_data %>%
    left_join(semi_age, by = ".draw") %>%  
    filter(REC_AGE > SEMI_AGE) %>%  
    group_by(.draw) %>%
    filter(DIG == max(DIG, na.rm = TRUE)) %>%  # Find max DIG proportion
    summarise(PEAK_DIG_AGE = first(REC_AGE), .groups = "drop")
  
  full_age <- call_data %>%
    left_join(peak_dig_age, by = ".draw") %>%  
    filter(REC_AGE > PEAK_DIG_AGE) %>%  
    filter(CC > DIG + PROP_DIFF) %>%  
    group_by(.draw) %>%
    summarise(FULL_AGE = min(REC_AGE, na.rm = TRUE), .groups = "drop")
  
  age_data <- semi_age %>%
    left_join(peak_dig_age, by = ".draw") %>%  # Merge peak DIG age
    left_join(full_age, by = ".draw")          # Merge full age
  
  # Ensure we have x samples per ID
  transition_ages <- age_data %>%
    slice_sample(n = num_samples, replace = F)  
  
  # Add metadata
  transition_ages <- transition_ages %>%
    mutate(
      ID = unique(ID_DATA$ID),
      TREATMENT = unique(ID_DATA$TREATMENT),
      SEX = unique(ID_DATA$SEX),
      WEIGHT = unique(ID_DATA$WEIGHT),
      COMP = unique(ID_DATA$COMP),
      GS = unique(ID_DATA$GS),
      RAIN = unique(ID_DATA$RAIN),
      MOTHER_ID = unique(ID_DATA$MOTHER_ID),
      LITTER_CODE = unique(ID_DATA$LITTER_CODE)
    )
  
  return(transition_ages)
}

# # ==============================================================================
# #  APPLY TO ALL INDIVIDUALS: Extract Posterior Call Proportions & Transitions
# # ==============================================================================
# # Create dataset with unique individuals but sample means for standardisation
# FLUT_ID_data <- PROP_data %>%
#   group_by(ID) %>%
#   summarise(
#     WEIGHT = mean(PROP_data$WEIGHT_DIFF_PER),
#     COMP = mean(PROP_data$COMP_NORM),
#     GS = mean(PROP_data$GROUPSIZE),
#     RAIN = mean(PROP_data$MonthlyRainfall),
#     WEIGHT_z = 0, COMP_z = 0, GS_z = 0, RAIN_z = 0,
#     MOTHER_ID = first(MOTHER_ID),
#     LITTER_CODE = first(LITTER_CODE),
#     ID = first(ID),
#     TREATMENT = as.factor(first(TREATMENT)),
#     SEX = first(SEX)
#   )
# 
# # Initialize an empty tibble to store results
# full_transition_data <- tibble()
# 
# # Loop through all subjects
# for (i in 1:nrow(FLUT_ID_data)) {
#   row <- FLUT_ID_data[i, ]
#   posterior_call_props <- extract_posterior_call_props(row)
#   transition_ages <- find_transition_ages(posterior_call_props, row, num_samples = NUM_SAMPLES_PER_ID)
#   
#   full_transition_data <- bind_rows(full_transition_data, transition_ages)
#   
#   # Print progress
#   percent_done <- round((i / nrow(FLUT_ID_data)) * 100, 2)
#   print(paste("Finished processing ID:", row$ID, "-", percent_done, "% complete"))
# }
# 
# 
# full_transition_data %>%
#   group_by(TREATMENT) %>%
#   summarise(n_samples = n()) %>%
#   arrange(desc(n_samples))
# 
# 
# # Save posterior transition ages
# saveRDS(full_transition_data, "MILESTONES_POSTERIOR.rds")
# 
# # Create sub datasets for each milestone
# # Load posterior transition age samples
# milestone_data <- readRDS("MILESTONES_POSTERIOR.rds")
# print(milestone_data %>%
#   group_by(TREATMENT, ID) %>%
#   summarise(n_samples = n()) %>%
#   arrange(desc(n_samples)), n=55)
# 
# 
# SEMI_data <- milestone_data %>%
#   select(ID, TREATMENT, SEX, WEIGHT, COMP, GS, RAIN, MOTHER_ID, LITTER_CODE, .draw, SEMI_AGE) %>%
#   drop_na(SEMI_AGE) %>%
#   rename(Milestone_Age = SEMI_AGE) %>%
#   mutate(Milestone_Type = "SEMI")
# saveRDS(SEMI_data, 'SEMI_data.rds')
# 
# PEAK_data <- milestone_data %>%
#   select(ID, TREATMENT, SEX, WEIGHT, COMP, GS, RAIN, MOTHER_ID, LITTER_CODE, .draw, PEAK_DIG_AGE) %>%
#   drop_na(PEAK_DIG_AGE) %>%
#   rename(Milestone_Age = PEAK_DIG_AGE) %>%
#   mutate(Milestone_Type = "PEAK")
# saveRDS(PEAK_data, 'PEAK_data.rds')
# 
# FULL_data <- milestone_data %>%
#   select(ID, TREATMENT, SEX, WEIGHT, COMP, GS, RAIN, MOTHER_ID, LITTER_CODE, .draw, FULL_AGE) %>%
#   drop_na(FULL_AGE) %>%
#   rename(Milestone_Age = FULL_AGE) %>%
#   mutate(Milestone_Type = "FULL")
# saveRDS(FULL_data, 'FULL_data.rds')


#  Transition from full to SEMI-dependence ####
SEMI_data <- readRDS('../data/SEMI_data.rds')

priors_SEMI <- c(set_prior("normal(0,1)", class = "Intercept"), set_prior("normal(0,1)", class='b')) 
SEMI_anova <- brm(formula = Milestone_Age  ~ (TREATMENT * SEX) + (1|ID),
                  data = SEMI_data,
                  family = lognormal(link='identity'), 
                  chains = 4, iter = 10000, warmup = 2500, seed = 42234223, control = list(max_treedepth = 15, adapt_delta=0.99),
                  save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                  prior = priors_SEMI, threads = threading(4),
                  file="SEMI_anova"
)


rm(priors_SEMI)

#### RESULTS: SEMI ####
SEMI_anova <- readRDS('../models/SEMI_anova.rds')

plot(SEMI_anova)
pp_check(SEMI_anova, ndraws=100)

loo_R2(SEMI_anova, moment_match=T)

bayes_R2(SEMI_anova)

performance::variance_decomposition(SEMI_anova)

summary(SEMI_anova)

describe_posterior(
  SEMI_anova,
  effects = "all", 
  component = "all",
  rope_range = rope_range(SEMI_anova),  
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)

# EMMs for TREATMENT:SEX
(treat_sex <- emmeans(SEMI_anova, pairwise ~ TREATMENT:SEX))

# in days:
exp_treat_sex <- as.data.frame(treat_sex$emmeans)
exp_treat_sex$emmean <- exp(exp_treat_sex$emmean)
exp_treat_sex$lower.HPD <- exp(exp_treat_sex$lower.HPD)
exp_treat_sex$upper.HPD <- exp(exp_treat_sex$upper.HPD)
print(exp_treat_sex)

pd(treat_sex)

p_significance(treat_sex, threshold = rope_range(SEMI_anova))

# Plot emmeans: TS ####
treat_sex <- emmeans(SEMI_anova, ~ TREATMENT:SEX)
emm_data <- as.data.frame(treat_sex)
emm_data$SEMI_age <- exp(emm_data$emmean)
emm_data$HPD_low <- exp(emm_data$lower.HPD)
emm_data$HPD_high <- exp(emm_data$upper.HPD)
emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("DC", "SC", 'DT'))

milestones_org <- SEMI_data %>%
  group_by(ID, TREATMENT, SEX) %>%
  summarise(SEMI_age = median(Milestone_Age, na.rm = TRUE), .groups = "drop")


# 700* 500: EMMs_SEMI_TS
ggplot(emm_data, aes(x = TREATMENT, y = SEMI_age, color = TREATMENT, shape = SEX)) +
  geom_point(data = milestones_org, aes(x = TREATMENT, y = SEMI_age, shape = SEX, color=TREATMENT), 
             position = position_dodge(width = 0.7), size = 3, alpha = 0.3) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.7), width = 0.5) +
  scale_color_okabe_ito(order= c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_y_continuous(n.breaks=10)+
  scale_x_discrete(labels= c('DC', 'SC', 'DT'))+
  labs(x = "Maternal status", y = "Age (days)\n") +
  theme_clean()

rm(treat_sex, milestones_org)

### PLOTS: SEMI ####
# coefficient plots ####
posterior_desc <- describe_posterior(
  SEMI_anova,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(SEMI_anova),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma
posterior_desc <- posterior_desc[-c(1, 7),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSC", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTDT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSC', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTDT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')

custom_order <- c( 'DT:Male','SC:Male',"Male", 
                   'DT','SC')

posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_SEMI 600*700
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1,3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()

#################
## DIG -> CC ####
FULL_data <- readRDS('../data/FULL_data.rds')

priors_FULL <- c(set_prior("normal(0,1)", class = "Intercept"), set_prior("normal(0,1)", class='b')) 
FULL_anova <- brm(formula = Milestone_Age ~ (TREATMENT * SEX) + (1|ID),
                  data = FULL_data,
                  family = lognormal(link='identity'), 
                  chains = 4, iter = 10000, warmup = 2500, seed = 42234223, control = list(max_treedepth = 15, adapt_delta=0.99),
                  save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                  prior = priors_FULL, threads = threading(4),
                  file="FULL_anova")
rm(priors_FULL)

### RESULTS: FULL ####
FULL_anova <- readRDS('../models/FULL_anova.rds')

plot(FULL_anova)
pp_check(FULL_anova, ndraws=100)

loo_R2(FULL_anova, moment_match=T)

bayes_R2(FULL_anova)

performance::variance_decomposition(FULL_anova)

summary(FULL_anova)

describe_posterior(
  FULL_anova,
  effects = "all", 
  component = "all",
  rope_range = rope_range(FULL_anova),  
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)

# EMMs for TREATMENT:SEX
(treat_sex <- emmeans(FULL_anova, pairwise ~ TREATMENT:SEX))

# in days:
exp_treat_sex <- as.data.frame(treat_sex$emmeans)
exp_treat_sex$emmean <- exp(exp_treat_sex$emmean)
exp_treat_sex$lower.HPD <- exp(exp_treat_sex$lower.HPD)
exp_treat_sex$upper.HPD <- exp(exp_treat_sex$upper.HPD)
print(exp_treat_sex)

pd(treat_sex)

p_significance(treat_sex, threshold = rope_range(FULL_anova))

# Plot emmeans: TS ####
treat_sex <- emmeans(FULL_anova, ~ TREATMENT:SEX)
emm_data <- as.data.frame(treat_sex)
emm_data$FULL_age <- exp(emm_data$emmean)
emm_data$HPD_low <- exp(emm_data$lower.HPD)
emm_data$HPD_high <- exp(emm_data$upper.HPD)
emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("DC", "SC", 'DT'))

milestones_org <- FULL_data %>%
  group_by(ID, TREATMENT, SEX) %>%
  summarise(FULL_age = median(Milestone_Age, na.rm = TRUE), .groups = "drop")

# 700* 500: EMMs_FULL_TS
ggplot(emm_data, aes(x = TREATMENT, y = FULL_age, color = TREATMENT, shape = SEX)) +
  geom_point(data = milestones_org, aes(x = TREATMENT, y = FULL_age, shape = SEX, color=TREATMENT), 
             position = position_dodge(width = 0.7), size = 3, alpha = 0.3) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.7), width = 0.5) +
  scale_color_okabe_ito(order= c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_y_continuous(n.breaks=10)+
  scale_x_discrete(labels= c('DC', 'SC', 'DT'))+
  labs(x = "Maternal status", y = "Age (days)\n") +
  theme_clean()

rm(treat_sex)

#### PLOTS: FULL ####
# coefficient plots #### 
posterior_desc <- describe_posterior(
  FULL_anova,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(FULL_anova),
  centrality = "median",
  dispersion = TRUE
)
# drop sigma
posterior_desc <- posterior_desc[-c(1,7),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSC", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTDT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTDT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSC', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')

custom_order <- c('DT:Male','SC:Male',"Male", 
                  'DT','SC')

posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_FULL 600*700
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  # ggtitle('Semi-dependence to full independence')+
  theme_clean()#+
# theme(legend.position="none")

rm(SEMI_data, FULL_data, SEMI_anova, FULL_anova)

#### PEAK DIG ####
PEAK_data <- readRDS('../data/PEAK_data.rds')

priors_DIG_peak <- c(set_prior("normal(0,1)", class = "Intercept"), set_prior("normal(0,1)", class='b')) 
DIG_peak_anova <- brm(formula = Milestone_Age ~ (TREATMENT * SEX) + (1|ID),
                      data = PEAK_data,
                      family = lognormal(link='identity'), 
                      chains = 4, iter = 10000, warmup = 2500, seed = 4223425, control = list(max_treedepth = 15, adapt_delta=0.9),
                      save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                      prior = priors_DIG_peak, threads = threading(4),
                      file="DIG_PEAK_anova"
)


rm(priors_DIG_peak)

### RESULTS PEAK ####
DIG_peak_anova <- readRDS('../models/DIG_PEAK_anova.rds')

plot(DIG_peak_anova)
pp_check(DIG_peak_anova, ndraws=100)

loo_R2(DIG_peak_anova, moment_match=T)

bayes_R2(DIG_peak_anova)


performance::variance_decomposition(DIG_peak_anova)

summary(DIG_peak_anova)
describe_posterior(
  DIG_peak_anova,
  effects = "all", 
  component = "all",
  rope_range = rope_range(DIG_peak_anova),  
  test = c("p_direction", "p_significance"),
  centrality = "all",
  dispersion = TRUE
)

# EMMs for TREATMENT:SEX
(treat_sex <- emmeans(DIG_peak_anova, pairwise ~ TREATMENT:SEX))

# in days:
exp_treat_sex <- as.data.frame(treat_sex$emmeans)
exp_treat_sex$emmean <- exp(exp_treat_sex$emmean)
exp_treat_sex$lower.HPD <- exp(exp_treat_sex$lower.HPD)
exp_treat_sex$upper.HPD <- exp(exp_treat_sex$upper.HPD)
print(exp_treat_sex)

pd(treat_sex)


p_significance(treat_sex, threshold = rope_range(DIG_peak_anova))

# Plot emmeans: TS ####
treat_sex <- emmeans(DIG_peak_anova, ~ TREATMENT:SEX)
emm_data <- as.data.frame(treat_sex)
emm_data$PEAK_age <- exp(emm_data$emmean)
emm_data$HPD_low <- exp(emm_data$lower.HPD)
emm_data$HPD_high <- exp(emm_data$upper.HPD)
emm_data$TREATMENT <- factor(emm_data$TREATMENT, levels = c("DC", "SC", 'DT'))

milestones_org <- PEAK_data %>%
  group_by(ID, TREATMENT, SEX) %>%
  summarise(PEAK_age = median(Milestone_Age, na.rm = TRUE), .groups = "drop")

# 700* 500: EMMs_PEAK_DIG_TS
ggplot(emm_data, aes(x = TREATMENT, y = PEAK_age, color = TREATMENT, shape = SEX)) +
  geom_point(data = milestones_org, aes(x = TREATMENT, y = PEAK_age, shape = SEX, color=TREATMENT), 
             position = position_dodge(width = 0.7), size = 3, alpha = 0.3) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.7), width = 0.5) +
  scale_color_okabe_ito(order= c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_y_continuous(n.breaks=10)+
  scale_x_discrete(labels= c('DC', 'SC', 'DT'))+
  labs(x = "Maternal status", y = "Age (days)\n") +
  theme_clean()

rm(treat_sex, milestones_org)

# PLOTS ####
# coefficient plots ####
posterior_desc <- describe_posterior(
  DIG_peak_anova,
  effects = "fixed",
  component = "all",
  rope_range = rope_range(DIG_peak_anova),
  centrality = "all",
  dispersion = TRUE
)
# drop sigma
posterior_desc <- posterior_desc[-c(1, 7),]

# clean up labels:
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'b_', '')
posterior_desc$TREATMENT <- ifelse(grepl("TREATMENTSC", posterior_desc$Parameter), "SC",
                                   ifelse(grepl("TREATMENTDT", posterior_desc$Parameter), "DT", "DC"))
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTSC', 'SC')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'TREATMENTDT', 'DT')
posterior_desc$Parameter <- str_replace_all(posterior_desc$Parameter, 'SEXM', 'Male')

custom_order <- c('DT:Male','SC:Male',"Male", 
                  'DT','SC')

posterior_desc$Parameter <- factor(posterior_desc$Parameter, levels = custom_order)
posterior_desc$TREATMENT <- factor(posterior_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_PEAK_DIG 600*700
ggplot(posterior_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()#+



#### ADD ON: PLOT SEMI; FULL; DIG PEAK ####
SEMI_anova <- readRDS('../models/SEMI_anova.rds')
FULL_anova <- readRDS('../models/FULL_anova.rds')
DIG_peak_anova <- readRDS('../models/DIG_PEAK_anova.rds')

SEMI_data <- readRDS('../data/SEMI_data.rds')
FULL_data <- readRDS('../data/FULL_data.rds')
PEAK_data <- readRDS('../data/PEAK_data.rds')

treat_sex <- emmeans(SEMI_anova, ~ TREATMENT:SEX)
emm_data_SEMI <- as.data.frame(treat_sex)
emm_data_SEMI$Age <- exp(emm_data_SEMI$emmean)
emm_data_SEMI$HPD_low <- exp(emm_data_SEMI$lower.HPD)
emm_data_SEMI$HPD_high <- exp(emm_data_SEMI$upper.HPD)
emm_data_SEMI$TREATMENT <- factor(emm_data_SEMI$TREATMENT, levels = c("DC", "SC", 'DT'))
emm_data_SEMI$Type <- 'SEMI'

treat_sex <- emmeans(FULL_anova, ~ TREATMENT:SEX)
emm_data_FULL <- as.data.frame(treat_sex)
emm_data_FULL$Age <- exp(emm_data_FULL$emmean)
emm_data_FULL$HPD_low <- exp(emm_data_FULL$lower.HPD)
emm_data_FULL$HPD_high <- exp(emm_data_FULL$upper.HPD)
emm_data_FULL$TREATMENT <- factor(emm_data_FULL$TREATMENT, levels = c("DC", "SC", 'DT'))
emm_data_FULL$Type <- 'FULL'

treat_sex <- emmeans(DIG_peak_anova, ~ TREATMENT:SEX)
emm_data_DIG <- as.data.frame(treat_sex)
emm_data_DIG$Age <- exp(emm_data_DIG$emmean)
emm_data_DIG$HPD_low <- exp(emm_data_DIG$lower.HPD)
emm_data_DIG$HPD_high <- exp(emm_data_DIG$upper.HPD)
emm_data_DIG$TREATMENT <- factor(emm_data_DIG$TREATMENT, levels = c("DC", "SC", 'DT'))
emm_data_DIG$Type <- 'DIG'

# if raw data should be plotted: SEMI, FULL and DIG_data needed!
# Combine raw data as well

semi_data_filtered <- SEMI_data
semi_data_filtered$Type <- "SEMI"

full_data_filtered <- FULL_data
full_data_filtered$Type <- "FULL"

dig_data_filtered <- PEAK_data
dig_data_filtered$Type <- "DIG"

# rename column in DIG data:
combined_raw_data <- rbind(semi_data_filtered, full_data_filtered, dig_data_filtered)
combined_raw_data$Type <- factor(combined_raw_data$Type, levels = c('SEMI', 'DIG', 'FULL'))

# calculate medians for raw data
combined_raw_data <- combined_raw_data %>%
  group_by(ID, TREATMENT, SEX, Type) %>%
  summarise(Milestone_Age = median(Milestone_Age, na.rm = TRUE), .groups = "drop")


data_combined <- rbind(emm_data_SEMI, emm_data_FULL, emm_data_DIG)
data_combined$Type <- as.factor(data_combined$Type)
data_combined$Type <- factor(data_combined$Type, levels = c('SEMI', 'DIG', 'FULL'))

facet_labels <- c('SEMI' = '(a)', 'DIG' = '(b)', 'FULL' = '(c)')

ggplot(data_combined, aes(x = TREATMENT, y = Age, color = TREATMENT, shape = SEX)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +  # Show estimates as points
  geom_errorbar(aes(ymin = HPD_low, ymax = HPD_high), position = position_dodge(width = 0.7), width = 0.5) +
  scale_color_okabe_ito(order= c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_shape_manual(name = "Sex", labels = c("female", "male"), values = c(16, 17)) + # Rename SEX to Sex and F/M to Female/Male
  scale_y_continuous(n.breaks=10)+
  scale_x_discrete(labels= c('DC', 'SC', 'DT'))+
  labs(x = "Maternal status", y = "Age (days)\n") +
  theme_clean()+
  facet_wrap(~ Type, nrow = 1, labeller = labeller(Type = facet_labels)) + # Wrap by transition
  # Filter first by Type, then take the first occurrence per ID
  geom_point(data = combined_raw_data %>% filter(Type == "SEMI") %>% group_by(ID) %>% slice(1) %>% ungroup(), 
             aes(x = TREATMENT, y = Milestone_Age, shape = SEX, color = TREATMENT), 
             position = position_dodge(width = 0.7), size = 3, alpha = 0.3) +
  geom_point(data = combined_raw_data %>% filter(Type == "DIG") %>% group_by(ID) %>% slice(1) %>% ungroup(), 
             aes(x = TREATMENT, y = Milestone_Age, shape = SEX, color = TREATMENT), 
             position = position_dodge(width = 0.7), size = 3, alpha = 0.3) +
  geom_point(data = combined_raw_data %>% filter(Type == "FULL") %>% group_by(ID) %>% slice(1) %>% ungroup(), 
             aes(x = TREATMENT, y = Milestone_Age, shape = SEX, color = TREATMENT), 
             position = position_dodge(width = 0.7), size = 3, alpha = 0.3) +
  theme(
    panel.spacing = unit(1.5, "lines")  # Adjust this value to increase or decrease spacing
  ) 

### cleanup ###
rm(milestone_data, B_prop, priors, transition_data, PROP_data)