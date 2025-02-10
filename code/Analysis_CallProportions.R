#### Analysis script for Flutamide study: CALL PROPORTIONS #####################
# Multivariate Bayesian Hierarchical Logistic Regression Model including
# Treatment (T), age (A), sex (S), Body mass offset (W),
# also checking for (normalised) competition score (C = H/P ratio & pups & adults),
# and (adult) group size 
#
# BWalkenhorst, 2024

#### SETUP ####
# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #includes hidden objects.
gc()

# load all necessary libraries
library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr) # annotate_figure
library(tidyverse)
library(tidybayes)
library(brms)
library(rstan)
library(bayestestR)
library(ggokabeito) # colour palette
library(emmeans) #
library(extrafont)
# use on first use
font_import()

set.seed(23)

PROP_data <- readRDS('../data/PROP_data.rds')

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

#### FUNCTIONS ####
get_org_value_for_z <- function(z_value, column){
  return ((z_value*sd(column)) + mean(column))   
}

get_z_value_for_org <- function(org_value, column){
  return ((org_value - mean(column))/sd(column))   
}

########################## BAYES MODEL ########################################
setwd('../models/')

priors <- c(
  # Intercepts
  set_prior("normal(0, 1)", class = "Intercept", resp = c("SumBEG", "SumDIG", "SumCC")),
  set_prior("normal(0, 1)", class = "b", resp = c("SumBEG", "SumDIG", "SumCC"))
)

# use TAG2, TAR, TAW, TACN2, TA, TAS
bf_REP <- bf(Sum_BEG | trials(Total_calls) ~ TREATMENT + AGE_z + SEX + WEIGHT_z + COMP_NORM_z + I(COMP_NORM_z^2) + GS_z + I(GS_z^2) + RAIN_z +
               TREATMENT:AGE_z + TREATMENT:SEX + TREATMENT:WEIGHT_z +TREATMENT:COMP_NORM_z + TREATMENT:I(COMP_NORM_z^2) + TREATMENT:GS_z +
               TREATMENT:I(GS_z^2) + TREATMENT:RAIN_z +
               AGE_z:SEX + AGE_z:WEIGHT_z + AGE_z:COMP_NORM_z + AGE_z:I(COMP_NORM_z^2) + AGE_z:GS_z + AGE_z:I(GS_z^2) + AGE_z:RAIN_z +
               TREATMENT:AGE_z:SEX + TREATMENT:AGE_z:WEIGHT_z + TREATMENT:AGE_z:COMP_NORM_z + TREATMENT:AGE_z:I(COMP_NORM_z^2) +
               TREATMENT:AGE_z:GS_z + TREATMENT:AGE_z:I(GS_z^2) + TREATMENT:AGE_z:RAIN_z + (1|LITTER_CODE/ID))
# use TAG, TACN2, TAS, TAR2, TA2, TAW
bf_DIG <-bf(Sum_DIG | trials(Total_calls) ~ TREATMENT + AGE_z + I(AGE_z^2) + SEX + WEIGHT_z + COMP_NORM_z + I(COMP_NORM_z^2) + GS_z + RAIN_z + I(RAIN_z^2) +
              TREATMENT:AGE_z + TREATMENT:I(AGE_z^2) + TREATMENT:SEX + TREATMENT:WEIGHT_z + TREATMENT:COMP_NORM_z + TREATMENT:I(COMP_NORM_z^2) +
              TREATMENT:GS_z + TREATMENT:RAIN_z + TREATMENT:I(RAIN_z^2) +
              AGE_z:SEX + AGE_z:WEIGHT_z + AGE_z:COMP_NORM_z + AGE_z:I(COMP_NORM_z^2) + AGE_z:GS_z + AGE_z:RAIN_z + AGE_z:I(RAIN_z^2) +
              I(AGE_z^2):SEX + I(AGE_z^2):WEIGHT_z + I(AGE_z^2):COMP_NORM_z + I(AGE_z^2):I(COMP_NORM_z^2) + I(AGE_z^2):GS_z + I(AGE_z^2):RAIN_z + I(AGE_z^2):I(RAIN_z^2) +
              TREATMENT:AGE_z:SEX + TREATMENT:AGE_z:WEIGHT_z + TREATMENT:AGE_z:COMP_NORM_z  + TREATMENT:AGE_z:I(COMP_NORM_z^2) + TREATMENT:AGE_z:GS_z +
              TREATMENT:AGE_z:RAIN_z + TREATMENT:AGE_z:I(RAIN_z^2) +
              TREATMENT:I(AGE_z^2):SEX + TREATMENT:I(AGE_z^2):WEIGHT_z + TREATMENT:I(AGE_z^2):COMP_NORM_z + TREATMENT:I(AGE_z^2):I(COMP_NORM_z^2) +
              TREATMENT:I(AGE_z^2):GS_z + TREATMENT:I(AGE_z^2):RAIN_z + TREATMENT:I(AGE_z^2):I(RAIN_z^2) + (1|LITTER_CODE/ID))
# use use TAR, TAG2, TACN2, TAS, TAW, TA
bf_CC <- bf(Sum_CC | trials(Total_calls) ~ TREATMENT + AGE_z + SEX + WEIGHT_z + COMP_NORM_z + I(COMP_NORM_z^2) + GS_z + I(GS_z^2)+ RAIN_z +
              TREATMENT:AGE_z + TREATMENT:SEX + TREATMENT:WEIGHT_z +TREATMENT:COMP_NORM_z + TREATMENT:I(COMP_NORM_z^2) + TREATMENT:GS_z +
              TREATMENT:I(GS_z^2) + TREATMENT:RAIN_z +
              AGE_z:SEX + AGE_z:WEIGHT_z + AGE_z:COMP_NORM_z + AGE_z:I(COMP_NORM_z^2) + AGE_z:GS_z + AGE_z:I(GS_z^2) + AGE_z:RAIN_z +
              TREATMENT:AGE_z:SEX + TREATMENT:AGE_z:WEIGHT_z + TREATMENT:AGE_z:COMP_NORM_z + TREATMENT:AGE_z:I(COMP_NORM_z^2) +
              TREATMENT:AGE_z:GS_z + TREATMENT:AGE_z:I(GS_z^2) + TREATMENT:AGE_z:RAIN_z + (1|LITTER_CODE/ID))

multivar_formula <- mvbrmsformula(bf_REP, bf_DIG, bf_CC)

B_prop <- brms::brm(formula = multivar_formula,
                    data = PROP_data, family = zero_inflated_beta_binomial(link='logit'),
                    chains = 4, iter = 5000, warmup = 1500, seed = 23, control = list(max_treedepth = 20, adapt_delta=0.99),
                    save_pars = save_pars(all=T), cores=4, backend = 'cmdstanr', init=0,
                    threads = threading(4),
                    prior = priors,
                    file="B_prop")

#### RESULTS: ####
B_prop <- readRDS('B_prop.rds')
options(max.print=999999)

summary(B_prop)

plot(B_prop)
pp_check(B_prop, ndraws = 100, resp='SumBEG')
pp_check(B_prop, ndraws = 100, resp='SumDIG')
pp_check(B_prop, ndraws = 100, resp='SumCC')

posterior <- describe_posterior(
  B_prop$fit,
  effects = "all", #fixed vs all ( vs random)
  component = "all",
  rope_range = c(-0.18, 0.18), 
  test = c("pd", "ps"),
  centrality = c('mean', 'median'),
  dispersion = TRUE#, 
 # parameters = 'SumBEG'
)
(output <- posterior[1:183,])


loo_R2(B_prop, moment_match=T)

bayes_R2(B_prop)


### Coefficient plots ####
posterior_desc <- output

REP_desc <- posterior_desc %>%
    filter(str_detect(Parameter, "SumBEG"))
DIG_desc <- posterior_desc %>%
  filter(str_detect(Parameter, "SumDIG"))
CC_desc <- posterior_desc %>%
  filter(str_detect(Parameter, "SumCC"))

# REP Coeff
REP_desc <- REP_desc[c(2:48),] 
# clean up labels:
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'b_SumBEG_', '')
REP_desc$TREATMENT <- ifelse(grepl("TREATMENTSC", REP_desc$Parameter), "SC",
                             ifelse(grepl("TREATMENTDT", REP_desc$Parameter), "DT", "DC"))
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'TREATMENTDT', 'DT')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'TREATMENTSC', 'SC')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'WEIGHT_z', 'Body mass offset')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'IGS_zE2', 'Group size^2')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'GS_z', 'Group size')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'AGE_z', 'Age')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'ICOMP_NORM_zE2', 'Competition^2')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'COMP_NORM_z', 'Competition')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'RAIN_z', 'Monthly rainfall')
REP_desc$Parameter <- str_replace_all(REP_desc$Parameter, 'SEXM', 'Male')

custom_order <- c('DT:Age:Monthly rainfall', 'SC:Age:Monthly rainfall','Age:Monthly rainfall', 
                  'DT:Age:Group size^2', 'SC:Age:Group size^2','Age:Group size^2', 
                  'DT:Age:Group size', 'SC:Age:Group size','Age:Group size', 
                  'DT:Age:Competition^2', 'SC:Age:Competition^2','Age:Competition^2',
                  'DT:Age:Competition', 'SC:Age:Competition','Age:Competition',
                  'DT:Age:Body mass offset', 'SC:Age:Body mass offset','Age:Body mass offset', 
                  'DT:Age:Male', 'SC:Age:Male','Age:Male', 
                  'DT:Monthly rainfall', 'SC:Monthly rainfall','Monthly rainfall', 
                  'DT:Group size^2', 'SC:Group size^2','Group size^2',
                  'DT:Group size', 'SC:Group size','Group size', 
                  'DT:Competition^2', 'SC:Competition^2','Competition^2',
                  'DT:Competition', 'SC:Competition','Competition', 
                  'DT:Body mass offset', 'SC:Body mass offset','Body mass offset',
                  'DT:Male', 'SC:Male','Male', 
                  'DT:Age', 'SC:Age', "Age",
                  'DT','SC')
# Update the order of TREATMENT factor levels
REP_desc$TREATMENT <- factor(REP_desc$TREATMENT, levels = c("DC", "SC", "DT"))
REP_desc$Parameter <- factor(REP_desc$Parameter, levels = custom_order)

# Coeff_REP 700*900
ggplot(REP_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()

# DIG coeff
DIG_desc <- DIG_desc[c(2:72),]

DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'b_SumDIG_', '')
DIG_desc$TREATMENT <- ifelse(grepl("TREATMENTSC", DIG_desc$Parameter), "SC",
                             ifelse(grepl("TREATMENTDT", DIG_desc$Parameter), "DT", "DC"))
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'TREATMENTDT', 'DT')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'TREATMENTSC', 'SC')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'WEIGHT_z', 'Body mass offset')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'IAGE_zE2', 'Age^2')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'IGS_zE2', 'Group size^2')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'GS_z', 'Group size')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'AGE_z', 'Age')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'ICOMP_NORM_zE2', 'Competition^2')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'COMP_NORM_z', 'Competition')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'IRAIN_zE2', 'Monthly rainfall^2')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'RAIN_z', 'Monthly rainfall')
DIG_desc$Parameter <- str_replace_all(DIG_desc$Parameter, 'SEXM', 'Male')

custom_order <- c('DT:Age^2:Monthly rainfall^2', 'SC:Age^2:Monthly rainfall^2','Age^2:Monthly rainfall^2', 
                  'DT:Age^2:Monthly rainfall', 'SC:Age^2:Monthly rainfall','Age^2:Monthly rainfall', 
                  'DT:Age^2:Group size^2', 'SC:Age^2:Group size^2','Age^2:Group size^2', 
                  'DT:Age^2:Group size', 'SC:Age^2:Group size','Age^2:Group size', 
                  'DT:Age^2:Competition^2', 'SC:Age^2:Competition^2','Age^2:Competition^2', 
                  'DT:Age^2:Competition', 'SC:Age^2:Competition','Age^2:Competition', 
                  'DT:Age^2:Body mass offset', 'SC:Age^2:Body mass offset','Age^2:Body mass offset', 
                  'DT:Age^2:Male', 'SC:Age^2:Male','Age^2:Male', 
                  'DT:Age:Monthly rainfall^2', 'SC:Age:Monthly rainfall^2','Age:Monthly rainfall^2', 
                  'DT:Age:Monthly rainfall', 'SC:Age:Monthly rainfall','Age:Monthly rainfall', 
                  'DT:Age:Group size^2', 'SC:Age:Group size^2','Age:Group size^2', 
                  'DT:Age:Group size', 'SC:Age:Group size','Age:Group size', 
                  'DT:Age:Competition^2', 'SC:Age:Competition^2','Age:Competition^2', 
                  'DT:Age:Competition', 'SC:Age:Competition','Age:Competition', 
                  'DT:Age:Body mass offset', 'SC:Age:Body mass offset','Age:Body mass offset', 
                  'DT:Age:Male', 'SC:Age:Male','Age:Male', 
                  'DT:Monthly rainfall^2', 'SC:Monthly rainfall^2','Monthly rainfall^2',
                  'DT:Monthly rainfall', 'SC:Monthly rainfall','Monthly rainfall',
                  'DT:Group size^2', 'SC:Group size^2', 'Group size^2', 
                  'DT:Group size', 'SC:Group size', 'Group size', 
                  'DT:Competition^2', 'SC:Competition^2','Competition^2',
                  'DT:Competition', 'SC:Competition','Competition', 
                  'DT:Body mass offset', 'SC:Body mass offset','Body mass offset', 
                  'DT:Male', 'SC:Male',"Male", 
                  'DT:Age^2', 'SC:Age^2',"Age^2", 
                  'DT:Age', 'SC:Age','Age', 
                  'DT', 'SC') 

DIG_desc$Parameter <- factor(DIG_desc$Parameter, levels = custom_order)
DIG_desc$TREATMENT <- factor(DIG_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_DIG 800*1100
ggplot(DIG_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()

# CC Coeff
CC_desc <- CC_desc[c(2:48),] 
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'b_SumCC_', '')
CC_desc$TREATMENT <- ifelse(grepl("TREATMENTSC", CC_desc$Parameter), "SC",
                             ifelse(grepl("TREATMENTDT", CC_desc$Parameter), "DT", "DC"))
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'TREATMENTDT', 'DT')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'TREATMENTSC', 'SC')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'WEIGHT_z', 'Body mass offset')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'IGS_zE2', 'Group size^2')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'GS_z', 'Group size')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'AGE_z', 'Age')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'ICOMP_NORM_zE2', 'Competition^2')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'COMP_NORM_z', 'Competition')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'RAIN_z', 'Monthly rainfall')
CC_desc$Parameter <- str_replace_all(CC_desc$Parameter, 'SEXM', 'Male')

custom_order <- c('DT:Age:Monthly rainfall', 'SC:Age:Monthly rainfall','Age:Monthly rainfall', 
                  'DT:Age:Group size^2', 'SC:Age:Group size^2','Age:Group size^2', 
                  'DT:Age:Group size', 'SC:Age:Group size','Age:Group size', 
                  'DT:Age:Competition^2', 'SC:Age:Competition^2','Age:Competition^2',
                  'DT:Age:Competition', 'SC:Age:Competition','Age:Competition',
                  'DT:Age:Body mass offset', 'SC:Age:Body mass offset','Age:Body mass offset', 
                  'DT:Age:Male', 'SC:Age:Male','Age:Male', 
                  'DT:Monthly rainfall', 'SC:Monthly rainfall','Monthly rainfall', 
                  'DT:Group size^2', 'SC:Group size^2','Group size^2',
                  'DT:Group size', 'SC:Group size','Group size', 
                  'DT:Competition^2', 'SC:Competition^2','Competition^2',
                  'DT:Competition', 'SC:Competition','Competition', 
                  'DT:Body mass offset', 'SC:Body mass offset','Body mass offset',
                  'DT:Male', 'SC:Male','Male', 
                  'DT:Age', 'SC:Age', "Age",
                  'DT','SC')

CC_desc$Parameter <- factor(CC_desc$Parameter, levels = custom_order)
CC_desc$TREATMENT <- factor(CC_desc$TREATMENT, levels = c("DC", "SC", "DT"))

# Coeff_CC 700*900
ggplot(CC_desc, aes(y = Parameter, x = Median, xmin = CI_low, xmax = CI_high, color=TREATMENT)) +
  geom_vline(xintercept = 0, color='grey', linetype = 'dotted', linewidth =1)+
  geom_point() +
  geom_errorbarh(height = 0) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT'))+
  labs(x = "Parameter value", y = "Parameters") +
  theme_clean()

rm(REP_desc, DIG_desc, CC_desc, output, posterior_desc, custom_order)

#### EMMs TREAT:SEX####
B_prop <- readRDS('../models/B_prop.rds')

### REP ####
# (REP_ROPE <- rope_range(B_REP_prop_TA)) # -0.1813799  0.1813799
REP_ROPE <- c(-0.18, 0.18)
(treat_sex <- emtrends(B_prop, pairwise ~ TREATMENT:SEX, var="AGE_z", resp='SumBEG'))

pd(treat_sex)

p_significance(treat_sex, threshold = REP_ROPE)


### DIG ####
# (DIG_ROPE <- rope_range(B_DIG_prop_TA2)) # -0.1813799  0.1813799
DIG_ROPE <- c(-0.18, 0.18)
(treat_sex <- emtrends(B_prop, pairwise ~ TREATMENT:SEX, var="AGE_z", resp='SumDIG', max.degree = 2))

pd(treat_sex)

p_significance(treat_sex, threshold = DIG_ROPE)

### CC ####
# (CC_ROPE <- rope_range(B_CC_prop_TA))# -0.1813799  0.1813799
CC_ROPE <- c(-0.18, 0.18)
(treat_sex <- emtrends(B_prop, pairwise ~ TREATMENT:SEX, var="AGE_z", resp='SumCC'))

pd(treat_sex)

p_significance(treat_sex, threshold = CC_ROPE)


rm(treat_sex, REP_ROPE, DIG_ROPE, CC_ROPE)

### EMMs at specific ages: 30, 75, 120 ###
ROPE <- c(-0.18, 0.18)
(REP_emms <- emmeans(B_prop,
                    pairwise ~ TREATMENT * SEX | AGE_z,
                    at = list(AGE_z = get_z_value_for_org(30, PROP_data$REC_AGE)),
                   resp='SumBEG'))
                  #  resp='SumBEG',type='response'))

p_significance(REP_emms, threshold = ROPE)

(DIG_emms <- emmeans(B_prop,
                     pairwise ~ TREATMENT * SEX | AGE_z,
                     at = list(AGE_z = get_z_value_for_org(75, PROP_data$REC_AGE)),
                      resp='SumDIG'))
#resp='SumDIG', type="response"))

pd(DIG_emms, rope_range = ROPE)

p_significance(DIG_emms, threshold = ROPE)

(CC_emms <- emmeans(B_prop,pairwise ~ TREATMENT * SEX | AGE_z,
                     at = list(AGE_z = get_z_value_for_org(120, PROP_data$REC_AGE)),
                     resp='SumCC'))
# resp='SumCC', type="response"))

pd(CC_emms, rope_range = ROPE)

p_significance(CC_emms, threshold = ROPE)

rm(ROPE, REP_emms, DIG_emms, CC_emms)

#### EMMs SOCIO-ECO####
B_prop <- readRDS('../models/B_prop.rds')

### REP ###
# # Body mass offset ###
(hyp_test_weight <- hypothesis(B_prop, c(
  "SumBEG_AGE_z:WEIGHT_z - SumBEG_TREATMENTSC:AGE_z:WEIGHT_z  < 0",   # DC vs. SC
  "SumBEG_AGE_z:WEIGHT_z  - SumBEG_TREATMENTDT:AGE_z:WEIGHT_z < 0",   # DC vs. DT
  "SumBEG_TREATMENTSC:AGE_z:WEIGHT_z - SumBEG_TREATMENTDT:AGE_z:WEIGHT_z < 0"  # SC vs. DT
)))

# Competition ###
(hyp_test_comp <- hypothesis(B_prop, c(
  "SumBEG_COMP_NORM_z - SumBEG_TREATMENTSC:COMP_NORM_z < 0",   # DC vs. SC
  "SumBEG_COMP_NORM_z - SumBEG_TREATMENTDT:COMP_NORM_z < 0",   # DC vs. DT
  "SumBEG_TREATMENTSC:COMP_NORM_z - SumBEG_TREATMENTDT:COMP_NORM_z > 0"  # SC vs. DT
)))

# GS ##
(hyp_test_gs <- hypothesis(B_prop, c(
  "SumBEG_AGE_z:GS_z - SumBEG_TREATMENTSC:AGE_z:GS_z < 0",   # DC vs. SC
  "SumBEG_AGE_z:GS_z - SumBEG_TREATMENTDT:AGE_z:GS_z < 0",   # DC vs. DT
  "SumBEG_TREATMENTSC:AGE_z:GS_z - SumBEG_TREATMENTDT:AGE_z:GS_z > 0"  # SC vs. DT
)))

# GS^2 ##
(hyp_test_gs2 <- hypothesis(B_prop, c(
  "SumBEG_AGE_z:IGS_zE2 - SumBEG_TREATMENTSC:AGE_z:IGS_zE2 < 0",   # DC vs. SC
  "SumBEG_AGE_z:IGS_zE2 - SumBEG_TREATMENTDT:AGE_z:IGS_zE2 < 0",   # DC vs. DT
  "SumBEG_TREATMENTSC:AGE_z:IGS_zE2 - SumBEG_TREATMENTDT:AGE_z:IGS_zE2 > 0"  # SC vs. DT
)))

# Rainfall ##
(hyp_test_rain <- hypothesis(B_prop, c(
  "SumBEG_RAIN_z - SumBEG_TREATMENTSC:RAIN_z < 0",   # DC vs. SC
  "SumBEG_RAIN_z - SumBEG_TREATMENTDT:RAIN_z < 0",   # DC vs. DT
  "SumBEG_TREATMENTSC:RAIN_z - SumBEG_TREATMENTDT:RAIN_z > 0"  # SC vs. DT
)))

### DIG ###
# Age:Competition2 ###
(hyp_test_comp2 <- hypothesis(B_prop, c(
  "SumDIG_AGE_z:ICOMP_NORM_zE2 - SumDIG_TREATMENTSC:AGE_z:ICOMP_NORM_zE2 > 0",   # DC vs. SC
  "SumDIG_AGE_z:ICOMP_NORM_zE2 - SumDIG_TREATMENTDT:AGE_z:ICOMP_NORM_zE2 > 0",   # DC vs. DT
  "SumDIG_TREATMENTSC:AGE_z:ICOMP_NORM_zE2 - SumDIG_TREATMENTDT:AGE_z:ICOMP_NORM_zE2 < 0"  # SC vs. DT
)))

# Age²:Competition ##
(hyp_test_comp <- hypothesis(B_prop, c(
  "SumDIG_IAGE_zE2:COMP_NORM_z - SumDIG_TREATMENTSC:IAGE_zE2:COMP_NORM_z > 0",   # DC vs. SC
  "SumDIG_IAGE_zE2:COMP_NORM_z - SumDIG_TREATMENTDT:IAGE_zE2:COMP_NORM_z > 0",   # DC vs. DT
  "SumDIG_TREATMENTSC:IAGE_zE2:COMP_NORM_z - SumDIG_TREATMENTDT:IAGE_zE2:COMP_NORM_z > 0"  # SC vs. DT
)))

# no age just comp
(hyp_test_comp <- hypothesis(B_prop, c(
  "SumDIG_COMP_NORM_z - SumDIG_TREATMENTSC:COMP_NORM_z < 0",   # DC vs. SC
  "SumDIG_COMP_NORM_z - SumDIG_TREATMENTDT:COMP_NORM_z < 0",   # DC vs. DT
  "SumDIG_TREATMENTSC:COMP_NORM_z - SumDIG_TREATMENTDT:COMP_NORM_z > 0"  # SC vs. DT
)))

# Age:Rainfall² ##
(hyp_test_rain2 <- hypothesis(B_prop, c(
  "SumDIG_AGE_z:IRAIN_zE2 - SumDIG_TREATMENTSC:AGE_z:IRAIN_zE2 > 0",   # DC vs. SC
  "SumDIG_AGE_z:IRAIN_zE2 - SumDIG_TREATMENTDT:AGE_z:IRAIN_zE2 > 0",   # DC vs. DT
  "SumDIG_TREATMENTSC:AGE_z:IRAIN_zE2 - SumDIG_TREATMENTDT:AGE_z:IRAIN_zE2 > 0"  # SC vs. DT
)))
 
         
# Age²:Rainfall ##
(hyp_test_rain <- hypothesis(B_prop, c(
  "SumDIG_IAGE_zE2:RAIN_z - SumDIG_TREATMENTSC:IAGE_zE2:RAIN_z < 0",   # DC vs. SC
  "SumDIG_IAGE_zE2:RAIN_z - SumDIG_TREATMENTDT:IAGE_zE2:RAIN_z > 0",   # DC vs. DT
  "SumDIG_TREATMENTSC:IAGE_zE2:RAIN_z - SumDIG_TREATMENTDT:IAGE_zE2:RAIN_z > 0"  # SC vs. DT
)))

### CC ###
# Age:Competition2 ###
(hyp_test_comp2 <- hypothesis(B_prop, c(
  "SumCC_AGE_z:ICOMP_NORM_zE2 - SumCC_TREATMENTSC:AGE_z:ICOMP_NORM_zE2 < 0",   # DC vs. SC
  "SumCC_AGE_z:ICOMP_NORM_zE2 - SumCC_TREATMENTDT:AGE_z:ICOMP_NORM_zE2 < 0",   # DC vs. DT
  "SumCC_TREATMENTSC:AGE_z:ICOMP_NORM_zE2 - SumCC_TREATMENTDT:AGE_z:ICOMP_NORM_zE2 < 0"  # SC vs. DT
)))


# GS² ##
(hyp_test_gs2 <- hypothesis(B_prop, c(
  "SumCC_AGE_z:IGS_zE2 - SumCC_TREATMENTSC:AGE_z:IGS_zE2 > 0",   # DC vs. SC
  "SumCC_AGE_z:IGS_zE2 - SumCC_TREATMENTDT:AGE_z:IGS_zE2 < 0",   # DC vs. DT
  "SumCC_TREATMENTSC:AGE_z:IGS_zE2 - SumCC_TREATMENTDT:AGE_z:IGS_zE2 < 0"  # SC vs. DT
)))


# Rainfall ##
(hyp_test_rain <- hypothesis(B_prop, c(
  "SumCC_AGE_z:RAIN_z - SumCC_TREATMENTSC:AGE_z:RAIN_z > 0",   # DC vs. SC
  "SumCC_AGE_z:RAIN_z - SumCC_TREATMENTDT:AGE_z:RAIN_z > 0",   # DC vs. DT
  "SumCC_TREATMENTSC:AGE_z:RAIN_z - SumCC_TREATMENTDT:AGE_z:RAIN_z > 0"  # SC vs. DT
)))

rm(hyp_test, hyp_test_comp, hyp_test_comp2, hyp_test_gs2, hyp_test_rain,
   hyp_test_weight)

#### MODEL PLOTS ####
B_prop <- readRDS('../models/B_prop.rds')
# get all needed values
sd_age <- sd(PROP_data$REC_AGE) #
mean_age <- mean(PROP_data$REC_AGE)

range(PROP_data$REC_AGE)# 31 130

rec_age_c <- seq(30, 130, by=1)
age_z_vals <- (rec_age_c - mean(PROP_data$REC_AGE))/sd(PROP_data$REC_AGE)
# 1 day steps not needed here 
age_z_vals <- seq(min(age_z_vals), max(age_z_vals), length.out = 20)

#predictions based on mean values
PROP_pred <- B_prop %>% 
  epred_draws(newdata = expand_grid(TREATMENT = levels(PROP_data$TREATMENT),
                                    SEX = levels(PROP_data$SEX),
                                    AGE_z = age_z_vals,
                                    WEIGHT_z = mean(PROP_data$WEIGHT_z),
                                    COMP_NORM_z = mean(PROP_data$COMP_NORM_z),
                                    GS_z = mean(PROP_data$GS_z),
                                    RAIN_z = mean(PROP_data$RAIN_z),
                                    Total_calls=1),
               re_formula = NA,   robust = T)

#unscale AGE_z values:
PROP_pred$REC_AGE <- PROP_pred$AGE_z * sd_age + mean_age
# ensure right format
PROP_pred$Call_prop <- PROP_pred$.epred 
PROP_pred$Call_type <- as.factor(PROP_pred$.category)
PROP_pred$SEX <- as.factor(PROP_pred$SEX)
PROP_pred$TREATMENT <- factor(PROP_pred$TREATMENT, levels = c("DC", "SC", "DT"))

ggplot(PROP_pred, aes(x = REC_AGE, y = Call_prop, color = TREATMENT, fill = TREATMENT, linetype=Call_type)) +  
  stat_lineribbon(.width = .95) + # shows uncertainty
  scale_linetype_manual(values = c("solid", "dotted", "twodash"), name = 'Call type', labels = c('REP', 'DIG', 'CC'),
                        guide = guide_legend(override.aes = list(color = "black"))) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Call proportion\n", 
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean()+
  facet_wrap(~SEX)


# just REP
REP_pred <- subset(PROP_pred, Call_type == 'SumBEG')
REP_pred$REP_prop <- REP_pred$Call_prop
#800*500: REP_TAS
ggplot(REP_pred, aes(x = REC_AGE, y = REP_prop, color = TREATMENT, fill = TREATMENT)) +  
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  geom_point(data = PROP_data, size = 1.5, alpha = 0.7) +  # raw data
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Repeat call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Repeat call proportion\n", 
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean()+
 facet_wrap(~SEX)

# just DIG
DIG_pred <- subset(PROP_pred, Call_type == 'SumDIG')
DIG_pred$DIG_prop <- DIG_pred$Call_prop
#800*500: DIG_TAS
ggplot(DIG_pred, aes(x = REC_AGE, y = DIG_prop, color = TREATMENT, fill = TREATMENT)) +  
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  geom_point(data = PROP_data, size = 1.5, alpha = 0.7) +  # raw data
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Digging call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Digging call proportion\n", 
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean()+
 # ggtitle('Sex')+
  facet_wrap(~SEX)

# just CC
CC_pred <- subset(PROP_pred, Call_type == 'SumCC')
CC_pred$CC_prop <- CC_pred$Call_prop
#800*500: CC_TAS
ggplot(CC_pred, aes(x = REC_AGE, y = CC_prop, color = TREATMENT, fill = TREATMENT)) +  
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  geom_point(data = PROP_data, size = 1.5, alpha = 0.7) +  # raw data
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Close call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Close call proportion\n", 
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean()+
  facet_wrap(~SEX)

# REP & DIG
REP_DIG_pred <- subset(PROP_pred, Call_type == 'SumBEG' | Call_type == 'SumDIG')
# 800*500: REP_DIG_TAS
ggplot(REP_DIG_pred, aes(x = REC_AGE, y = Call_prop, color = TREATMENT, fill = TREATMENT, linetype=Call_type)) +  
  stat_lineribbon(.width = .95) +
  scale_linetype_manual(values = c("solid", "twodash"), name = 'Call type', labels = c('REP', 'DIG'),
                        guide = guide_legend(override.aes = list(color = "black"))) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Call proportion\n", 
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean()+
  facet_wrap(~SEX)

# DIG & CC
DIG_CC_pred <- subset(PROP_pred, Call_type == 'SumDIG' | Call_type == 'SumCC')
# 800*500: DIG_CC_TAS
ggplot(DIG_CC_pred, aes(x = REC_AGE, y = Call_prop, color = TREATMENT, fill = TREATMENT, linetype=Call_type)) +  
  stat_lineribbon(.width = .95) +
  scale_linetype_manual(values = c("solid", "twodash"), name = 'Call type', labels = c('DIG', 'CC'),
                        guide = guide_legend(override.aes = list(color = "black"))) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Call proportion\n", 
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean()+
  facet_wrap(~SEX)

rm(REP_pred, DIG_pred, CC_pred, REP_DIG_pred, DIG_CC_pred)

### Normalised figures ####
# normalise predictions to ensure they sum up to 1
PROP_pred_sum <- PROP_pred %>%
  group_by(REC_AGE, TREATMENT, SEX, .draw) %>%
  mutate(Sum_prop = sum(Call_prop))

PROP_pred_normalized <- PROP_pred_sum %>%
  mutate(Normalized_prop = Call_prop / Sum_prop)

rm(PROP_pred_sum, PROP_pred)

#900*500: TAS (_re/_na)
ggplot(PROP_pred_normalized, aes(x = REC_AGE, y = Normalized_prop, color = TREATMENT, fill = TREATMENT, linetype=Call_type)) +  
  stat_lineribbon(.width = .95) + # shows uncertainty
  scale_linetype_manual(values = c("solid", "dotted", "twodash"), name = 'Call type', labels = c('REP', 'DIG', 'CC'),
                        guide = guide_legend(override.aes = list(color = "black"))) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Call proportion\n", 
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean()+
  facet_wrap(~SEX)

# just REP
REP_pred <- subset(PROP_pred_normalized, Call_type == 'SumBEG')
REP_pred$REP_prop <- REP_pred$Normalized_prop
#700*500: REP_TAS
ggplot(REP_pred, aes(x = REC_AGE, y = REP_prop, color = TREATMENT, fill = TREATMENT)) +  
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  geom_point(data = PROP_data, size = 1.5, alpha = 0.7) +  # raw data
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Repeat call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Repeat call proportion\n", 
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean()+
  facet_wrap(~SEX)

# just DIG
DIG_pred <- subset(PROP_pred_normalized, Call_type == 'SumDIG')
DIG_pred$DIG_prop <- DIG_pred$Normalized_prop
#700*500: DIG_TAS
ggplot(DIG_pred, aes(x = REC_AGE, y = DIG_prop, color = TREATMENT, fill = TREATMENT)) +  
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  geom_point(data = PROP_data, size = 1.5, alpha = 0.7) +  # raw data
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Digging call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Digging call proportion\n", 
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean()+
  facet_wrap(~SEX)

# just CC
CC_pred <- subset(PROP_pred_normalized, Call_type == 'SumCC')
CC_pred$CC_prop <- CC_pred$Normalized_prop
#700*500: CC_TAS
ggplot(CC_pred, aes(x = REC_AGE, y = CC_prop, color = TREATMENT, fill = TREATMENT)) +  
  stat_lineribbon(.width = .95) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  geom_point(data = PROP_data, size = 1.5, alpha = 0.7) +  # raw data
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Close call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Close call proportion\n", 
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean()+
  facet_wrap(~SEX)

# REP & DIG
REP_DIG_pred <- subset(PROP_pred_normalized, Call_type == 'SumBEG' | Call_type == 'SumDIG')
# 800*500: REP_DIG_TAS
ggplot(REP_DIG_pred, aes(x = REC_AGE, y = Normalized_prop, color = TREATMENT, fill = TREATMENT, linetype=Call_type)) +  
  stat_lineribbon(.width = .95) +
  scale_linetype_manual(values = c("solid", "twodash"), name = 'Call type', labels = c('REP', 'DIG'),
                        guide = guide_legend(override.aes = list(color = "black"))) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Call proportion\n", 
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean()+
  facet_wrap(~SEX)

# DIG & CC
DIG_CC_pred <- subset(PROP_pred_normalized, Call_type == 'SumDIG' | Call_type == 'SumCC')
# 800*500: DIG_CC_TAS
ggplot(DIG_CC_pred, aes(x = REC_AGE, y = Normalized_prop, color = TREATMENT, fill = TREATMENT, linetype=Call_type)) +  
  stat_lineribbon(.width = .95) +
  scale_linetype_manual(values = c("solid", "twodash"), name = 'Call type', labels = c('DIG', 'CC'),
                        guide = guide_legend(override.aes = list(color = "black"))) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(2, 1, 3), alpha = 0.3, name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  labs(x = "Age (days)", y = "Call proportion\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Call proportion\n", 
                     limits = c(0.0, 1.0),
                     breaks = c(0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) +
  theme_clean()+
  facet_wrap(~SEX)

rm(REP_pred, DIG_pred, CC_pred, REP_DIG_pred, DIG_CC_pred, summary_stats,
   age_z_vals, rec_age_c, mean_age, sd_age, PROP_pred_normalized)


