###### Model weights of FLUT  litters and determine offset to expected weight#####
#### BWalkenhorst 2024 ####

# clear everything and use garbage collector at start
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() 

# Load required libraries
library(brms)
library(readxl)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggokabeito) # colour palette
library(zoo) # rollapply
library(writexl)
library(bayestestR)
library(performance)
library(tidybayes)
library(extrafont)
# use on first use only
#font_import()

# set seed to duplicate results
set.seed(42)
# half normal prior and weak normal for intercept
priors_halfnormal <- c(set_prior('normal(0,0.5)', class = 'b', lb = 0), set_prior("normal(0,1)", class = "Intercept"))

# #Load and combine data ONLY NEEDED ONCE
# # weather data
# weather_data <- read_excel(
#   "sheets/weather/WEATHER20112014.xlsx"
# )
# 
# # subject weights
# weight_data <- read_excel(
#   "sheets/bodymass/FLUT_weights.xlsx"
# )
# # study data
# flut_data <- read_excel(
#   "sheets/4_FLUT_OFFSPRING_BW.xlsx",
# )
# # Data preprocessing
# flut_data_in <- flut_data %>%
#   select(File, REC_DATE, ID, TREATMENT, REC_AGE_D, DOB, SEX) %>%
#   mutate(REC_DATE = as.POSIXct(REC_DATE), TREATMENT = as.factor(TREATMENT),  DOB =as.POSIXct(DOB), SEX = as.factor(SEX))
# 
# weight_data$WeightTime <- as.POSIXct(weight_data$WeightTime)
# weight_data$ID <- weight_data$IndividID
# 
# # Filter to retain only the earliest weight per day for each ID
# earliest_weights <- weight_data %>%
#   group_by(ID, WeightDate) %>%                       # Group by ID and WeightDate
#   filter(WeightTime == min(WeightTime)) %>%          # Keep the row with the earliest WeightTime
#   ungroup()                                          # Ungroup to return to the original structure
# 
# # remove duplicate entries
# earliest_weights <- distinct(earliest_weights)
# 
# first_flut_data <- flut_data_in %>%
#   group_by(ID) %>%
#   slice(1) %>%
#   ungroup()
# 
# joined_data <- inner_join(earliest_weights, first_flut_data %>% select(ID, TREATMENT, REC_AGE_D, DOB, SEX), by = "ID")%>%
#   filter(AGE_DOB <= 140) # as max 136 needed for later avg calculations
# 
# table(joined_data$SEX)
# # F    M 
# # 1514 1539 
# 
# #sample_df <- joined_data[!duplicated(joined_data$ID),]
# # # 55 individuals
# #table(sample_df$SEX)
# # # F  M
# # # 29 26
# 
# # CUMULATIVE RAINFALL !!!
# weather_data <- weather_data %>%
#   arrange(Timestamp) %>%
#   mutate(Rainfall30D = rollapply(Rainfall24, width = 30, FUN = sum, align = "right", fill = NA))
# 
# weather_data <- weather_data %>%
#   ungroup()
# 
# # add rainfall to weight dates
# final_data <- joined_data %>%
#   left_join(weather_data %>% select(Timestamp, Rainfall30D), by = c("WeightDate" = "Timestamp"))

#saveRDS(final_data, '..//data/BM_data.rds')

final_data <- readRDS('../data/BM_data.rds')

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

#### MODEL ####
B_BODY <- brms::brm(formula = scale(Weight) ~ scale(AGE_DOB) + SEX + scale(Rainfall30D) + (1|ID), 
                        data = final_data, family = gaussian(link='identity'),
                        chains = 4, iter = 5000, warmup = 1500, cores = 4, backend = "cmdstanr", 
                        prior = priors_halfnormal,
                        seed = 23542235,
                        # save_pars = save_pars(all = TRUE),
                        control = list(max_treedepth = 15, adapt_delta=0.999), 
                        init=0, 
                        threads = threading(4),
                        file ="BodyMass.rds")

#### Model details ####
summary(B_BODY)

plot(B_BODY)
pp_check(B_BODY, ndraws=100)

# get the rope range  = -0.1 * SDy, 0.1 * SDy
# as its scaled, sd = 1 so -0.1 and 0.1 it is!
ropeRange <- c(-0.1* sd(scale(final_data$Weight)), 0.1 * sd(scale(final_data$Weight)))

describe_posterior(
  B_BODY,
  effects = "all", #fixed vs all (for random effects)
  component = "all",
  rope_range = ropeRange,  
  test = c("p_direction", "p_significance", "rope"),
  centrality = "all",
  dispersion = TRUE
)

loo_R2(B_BODY)


variance_decomposition(B_BODY)

### predict weight for FLUT pups ####
B_BODY <- readRDS("B_BODY.rds")

# # FLUT data with average weight and rainfall
# flut_data <- read_excel(
#   "G:/My Drive/Uni/UZH/Projects/KMP_PUPS/1stGen/DEV/sheets/5_FLUT_OFFSPRING_BW_Weight_Rain.xlsx",
# )
# flut_data <- flut_data %>%
#   mutate(TREATMENT = as.factor(TREATMENT), SEX = as.factor(SEX), ID = as.factor(ID))
# 
# # save flut data
# saveRDS(flut_data, 'BodyMass_FLUT_data.rds')
# 

flut_data <- readRDS('..//data/BodyMass_FLUT_data.rds')

copy_flut <- flut_data
copy_flut$AGE_DOB <- copy_flut$REC_AGE_D
copy_flut$Rainfall30D <- copy_flut$MonthlyRainfall

# get the predictions from the model
weight_pred <- predict(WEIGHT_AGE_RAIN, newdata = copy_flut, seed=23, allow_new_levels=F)

#saveRDS(weight_pred, "FLUT_BodyMass_pred.rds")
rm(copy_flut)

#load predictions if not calculated
weight_pred <- readRDS("..//data/FLUT_BodyMass_pred.rds")

flut_data$WEIGHT_PRED <- as.numeric(weight_pred[,1])* sd(flut_data$AvgWeight, na.rm = T) + mean(flut_data$AvgWeight, na.rm = T)
flut_data$WEIGHT_CI.lo <- as.numeric(weight_pred[,3])* sd(flut_data$AvgWeight, na.rm = T) + mean(flut_data$AvgWeight, na.rm = T)
flut_data$WEIGHT_CI.hi <- as.numeric(weight_pred[,4])* sd(flut_data$AvgWeight, na.rm = T) + mean(flut_data$AvgWeight, na.rm = T)

#plot predictions and actual average weight (jitter)
ggplot(flut_data, aes(x = REC_AGE_D, y = WEIGHT_PRED, color = TREATMENT)) +
  geom_jitter(alpha = 1, aes(shape = "predicted")) +
  geom_smooth(method = 'glm',  alpha = 0.3) +
  geom_jitter(aes(y = AvgWeight, shape = "average")) +
  geom_segment(aes(x = REC_AGE_D, xend = REC_AGE_D, y = AvgWeight, yend = WEIGHT_PRED), linetype = "dashed") +
  scale_fill_okabe_ito() +
  scale_color_okabe_ito(order = c(2, 3, 1), name = "Maternal\nstatus", labels = c('DC', 'DT', 'SC')) +
  labs(x = "Age (days)", y = "Body mass (g)") +
  scale_x_continuous(breaks = seq(0, 130, 10)) +
  scale_y_continuous(n.breaks = 10) +
  theme_clean() +
  guides(fill = FALSE, color = guide_legend(title = "Maternal\nstatus"),
         shape = guide_legend(title = "Weight"))

# epred
min(final_data$AGE_DOB)
max(final_data$AGE_DOB)

# expectation of the posterior predictive distribution: epred
prediction130D_data <- WEIGHT_AGE_RAIN %>% 
  epred_draws(newdata = expand_grid(AGE_DOB=seq(min(final_data$AGE_DOB), max(final_data$AGE_DOB), by=5),
                                    SEX = levels(final_data$SEX),
                                    Rainfall30D = mean(final_data$Rainfall30D)),
              re_formula = NULL, allow_new_levels=T, robust=T)

prediction130D_data$Weight <- prediction130D_data$.epred * sd(final_data$Weight) + mean(final_data$Weight)
flut_data$TREATMENT <- factor(flut_data$TREATMENT, levels = c("CTRL", "SUB", "FLUT"))

ggplot(prediction130D_data, aes(x = AGE_DOB, y = Weight)) +  
  stat_lineribbon(.width = .95, alpha = 0.5) +
  geom_segment(data = flut_data, aes(x = REC_AGE_D, xend = REC_AGE_D, y = AvgWeight, yend = WEIGHT_PRED, color = TREATMENT), linetype = "dashed", alpha = 0.7) +
  geom_jitter(data = flut_data, aes(x = REC_AGE_D, y = AvgWeight, colour = TREATMENT, shape = 'average'), alpha = 1, size = 1.5) +
  geom_jitter(data = flut_data, aes(x = REC_AGE_D, y = WEIGHT_PRED, colour = TREATMENT, shape = 'predicted'), alpha = 1, size = 2) +
  scale_color_okabe_ito(order = c(2, 1, 3), name = "Maternal\nstatus", labels = c('DC', 'SC', 'DT')) +
  scale_fill_okabe_ito(order = c(8), alpha = 0.2, guide = 'none') +
  labs(x = "Age (days)", y = "Body mass (g)\n") +
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_shape_manual(name = "Body mass", values = c('average' = 16, 'predicted' = 3)) +
  theme_clean()+
  guides(fill=FALSE)

############## determine OFFSET ####
flut_data$WEIGHT_DIFF <- flut_data$AvgWeight - flut_data$WEIGHT_PRED 

#use percentage
flut_data$WEIGHT_DIFF_PER <-as.numeric(((flut_data$AvgWeight - flut_data$WEIGHT_PRED)/flut_data$WEIGHT_PRED) *100)

# Plot offset
ggplot(flut_data, aes(x = REC_AGE_D, y = WEIGHT_DIFF_PER, color = TREATMENT)) +
  geom_point() +
 geom_smooth(method='glm', )+
  scale_x_continuous(breaks = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130), guide = guide_axis(angle = 45)) +
  scale_color_okabe_ito(order = c(2, 3, 1), name = "Maternal\nstatus", labels = c('DC', 'DT', 'SC')) +
  labs(x = "Age (days)",
       y = "Body mass offset (%)",
       color = "Treatment") +
  theme_clean()

# save the weight offset as excel sheet if needed
write_xlsx(flut_data, "/data/Filenam.xlsx")



