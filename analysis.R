rm(list = ls())

# ----------------------------------------------- packages -------------------------------------------------- #

# needed libraries

packages_needed <- c("readxl", "writexl", "dplyr","ggplot2", "corrplot", "tidyverse", 
                     "tidyquant", "ggdist", "ggthemes", "rstatix", "lme4", "fitdistrplus",
                     "performance", "gridExtra", "simr", "pastecs", "flexplot", 
                     "influence.ME", "MuMIn", "gridGraphics", "DHARMa", "lmerTest", "sjPlot", 
                     "ggeffects", "TOSTER", "broom.mixed", "emmeans", "purrr", "readr", 
                     "car", "lmtest", "irr", "MASS")


lapply(packages_needed, FUN = require, character.only = T)

get_package_info <- function(package_name) {
  if (requireNamespace(package_name, quietly = TRUE)) {
    pkg_desc <- packageDescription(package_name)
    
    package_info <- data.frame(
      Package = package_name,
      Description = ifelse("Title" %in% names(pkg_desc), pkg_desc[["Title"]], NA),
      Version = ifelse("Version" %in% names(pkg_desc), pkg_desc[["Version"]], NA)
    )
    return(package_info)
  } else {
    return(NULL)
  }
}

# Get information for each loaded package
package_info <- lapply(packages_needed, get_package_info)

# Filter out NULL values (packages that weren't loaded)
package_info <- Filter(Negate(is.null), package_info)

# Combine the information into a single data frame
package_info_df <- do.call(rbind, package_info)

# Make a table
writexl::write_xlsx(package_info_df, "PAPE_isometric_package_info.xlsx")
packages <- read_excel("PAPE_isometric_package_info.xlsx")

packages %>% as_tibble() %>% print(n = nrow(packages))

# --------------------------------------------------------------------------------------------------------- #
# --------------------------- setting the set.seed to returns the same results each run ------------------- #
# --------------------------------------------------------------------------------------------------------- #

set.seed(1234)

# --------------------------------------------------------------------------------------------------------- #
# ---------------------------------------- Randomization of testing sessions ------------------------------ #
# --------------------------------------------------------------------------------------------------------- #


participants_1 <- data.frame(ID = 1:30)

# Define a function to generate a random permutation of integers from 1 to 2
generate_permutation <- function() {
  sample(1:2)
}

# Define a vector of strings corresponding to the integers 1 to 2
replacement_strings_state <- c("baseline", "PAPE")
replacement_strings_1RM <- c("BP", "SMP")

# Apply the function to generate two new columns
participants_1[, 2:3] <- t(replicate(nrow(participants_1), generate_permutation()))

participants_1 <- participants_1 %>%
  rename(testing_1_1RM = V2, testing_2_1RM = V3)

# print(participants_1)

# Replace integers with strings in the new columns
for (i in 2:3) {
  participants_1[[i]] <- replacement_strings_1RM[participants_1[[i]]]
}

# Apply the function to generate two new columns
participants_1[, 4:5] <- t(replicate(nrow(participants_1), generate_permutation()))

participants_1 <- participants_1 %>%
  rename(testing_3_BP = V4, testing_4_BP = V5)

# print(participants_1)

# Replace integers with strings in the new columns
for (i in 4:5) {
  participants_1[[i]] <- replacement_strings_state[participants_1[[i]]]
}

# Print the resulting data frame
# print(participants_1)

# Apply the function to generate two new columns
participants_1[, 6:7] <- t(replicate(nrow(participants_1), generate_permutation()))

participants_1 <- participants_1 %>%
  rename(testing_5_SMP = V6, testing_6_SMP = V7)

for (i in 6:7) {
  participants_1[[i]] <- replacement_strings_state[participants_1[[i]]]
}

# Print the resulting data frame
print(participants_1)

setwd("E:/data/Statistics/Data")

write_xlsx(participants_1, "participants_random_select.xlsx")



# ---------------------------------------------------------------------------------------------------------- #
# --------------------------------------------- Data simulation -------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------- #

# IBP #
# We need to set a 1RM,
# NO-PAPE (2/5/7/10 time points),
# PAPE (2/5/7/10 time points)

# ISMP #
# We need to set a 1RM,
# NO-PAPE (2/5/7/10 time points),
# PAPE (2/5/7/10 time points)


# ----------------------------------- Isometric bench press ------------------------------------------------ #

# Based on the study of Kilduff et al. (2002) 'Effects of creatine on isometric bench-press performance in resistance-trained humans'
# SD = 255 N
# set fixed effect parameters
IBP_PFO_beta_0 <- 815 # intercept; i.e., the grand mean of PFO
# setting the slope parameter based on the Cohens d 0.3
# (M2 = 0.3 * 255N = 76.5N)...891.5-815 = 76.5
IBP_PFO_beta_1 <- 76.5 # slope; i.e, effect of category BASE/PAPE 

# set random effect parameters
IBP_PFO_tau_0 <- 75 # by-subject random intercept sd
IBP_PFO_omega_0 <- 10 # by-item random intercept sd (expect that BASE and PAPE have both the same sd)

# set more random effect and error parameters
IBP_PFO_tau_1 <- 30 # by-subject random slope sd
IBP_PFO_rho <- 0.2 # correlation between intercept and slope
IBP_PFO_sigma <- 20 # residual (error) sd (twice the size of the by-subject random intercept SD)

# simulating the sampling process

# set number of subjects and items
IBP_PFO_n_subj <- 30 # number of subjects
IBP_PFO_n_BASE <- 12 # number of (4 time points for BASE) 
IBP_PFO_n_PAPE <- 12 # number of (4 time points for PAPE) 


IBP_PFO_items <- data.frame(
  IBP_PFO_item_id = seq_len(IBP_PFO_n_BASE + IBP_PFO_n_PAPE),
  IBP_PFO_state = rep(c("BASE", "PAPE"), c(IBP_PFO_n_BASE, IBP_PFO_n_PAPE)),
  IBP_PFO_O_0i = rnorm(IBP_PFO_n_BASE + IBP_PFO_n_PAPE, mean = 0, sd = IBP_PFO_omega_0)
)

# effect-code category
IBP_PFO_items$IBP_PFO_X_i <- dplyr::recode(IBP_PFO_items$IBP_PFO_state, "BASE" = -0.5, "PAPE" = +0.5)


# simulate a sample of subjects
# calculate random intercept / random slope covariance
IBP_PFO_covar <- IBP_PFO_rho * IBP_PFO_tau_0 * IBP_PFO_tau_1

# put values into variance-covariance matrix
IBP_PFO_cov_mx  <- matrix(
  c(IBP_PFO_tau_0^2, IBP_PFO_covar,
    IBP_PFO_covar, IBP_PFO_tau_1^2),
  nrow = 2, byrow = TRUE)

# generate the by-subject random effects
IBP_PFO_subject_rfx <- MASS::mvrnorm(n = IBP_PFO_n_subj,
                                 mu = c(IBP_PFO_T_0s = 0, IBP_PFO_T_1s = 0),
                                 Sigma = IBP_PFO_cov_mx)

# combine with subject IDs
IBP_PFO_subjects <- data.frame(IBP_PFO_subj_id = seq_len(IBP_PFO_n_subj),
                               IBP_PFO_subject_rfx)

# simulate a sample of subjects
# sample from a multivariate random distribution 
IBP_PFO_subjects <- faux::rnorm_multi(
  n = IBP_PFO_n_subj, 
  mu = 0, # means for random effects are always 0
  sd = c(IBP_PFO_tau_0, IBP_PFO_tau_1), # set SDs
  r = IBP_PFO_rho, # set correlation, see ?faux::rnorm_multi
  varnames = c("IBP_PFO_T_0s", "IBP_PFO_T_1s")
)

# add subject IDs
IBP_PFO_subjects$IBP_PFO_subj_id <- seq_len(IBP_PFO_n_subj)


# cross subject and item IDs; add an error term
# nrow(.) is the number of rows in the table
IBP_PFO_trials <- crossing(IBP_PFO_subjects, IBP_PFO_items)  %>%
  mutate(IBP_PFO_e_si = rnorm(nrow(.), mean = 0, sd = IBP_PFO_sigma)) %>%
  dplyr::select(IBP_PFO_subj_id, IBP_PFO_item_id, IBP_PFO_state, IBP_PFO_X_i, everything())

IBP_PFO_trials <- IBP_PFO_trials[order(IBP_PFO_trials$IBP_PFO_subj_id, decreasing = F),]


# calculate the response variable
dat_sim_IBP_PFO <- IBP_PFO_trials %>%
  mutate(IBP_PFO = IBP_PFO_beta_0 + IBP_PFO_T_0s + IBP_PFO_O_0i + (IBP_PFO_beta_1 + IBP_PFO_T_1s) * IBP_PFO_X_i + IBP_PFO_e_si) %>%
  dplyr::select(IBP_PFO_subj_id, IBP_PFO_item_id, IBP_PFO_state, IBP_PFO_X_i, IBP_PFO)

View(dat_sim_IBP_PFO)


# -------------------------------- Isometric seated military press ----------------------------------------- #

# Based on the study of Vera et al. (2020) 'Intraocular Pressure Responses to Four Different Isometric Exercises in Men and Women'
#SD = 33.32 N
# set fixed effect parameters
ISMP_PFO_beta_0 <- 146 # intercept; i.e., the grand mean of PFO
# setting the slope parameter based on the Cohens d 0.3
# (M2 = 0.3 * 33.32N = 10.0N)...156-146 = 10.0N
ISMP_PFO_beta_1 <- 10 # slope; i.e, effect of category BASE/PAPE 

# set random effect parameters
ISMP_PFO_tau_0 <- 20 # by-subject random intercept sd
ISMP_PFO_omega_0 <- 0 # by-item random intercept sd (expect that BASE and PAPE have both the same sd)

# set more random effect and error parameters
ISMP_PFO_tau_1 <- 7 # by-subject random slope sd
ISMP_PFO_rho <- 0.2 # correlation between intercept and slope
ISMP_PFO_sigma <- 5 # residual (error) sd (twice the size of the by-subject random intercept SD)

# simulating the sampling process

# set number of subjects and items
ISMP_PFO_n_subj <- 30 # number of subjects
ISMP_PFO_n_BASE <- 12 # number of (3*4 BASE) 
ISMP_PFO_n_PAPE <- 12 # number of (3*4 PAPE) 


# total number of items = BASE + PAPE
ISMP_PFO_items <- data.frame(
  ISMP_PFO_item_id = seq_len(ISMP_PFO_n_BASE + ISMP_PFO_n_PAPE),
  ISMP_PFO_state = rep(c("BASE", "PAPE"), c(ISMP_PFO_n_BASE, ISMP_PFO_n_PAPE)),
  ISMP_PFO_O_0i = rnorm(ISMP_PFO_n_BASE + ISMP_PFO_n_PAPE, mean = 0, sd = ISMP_PFO_omega_0)
)

# effect-code category
ISMP_PFO_items$ISMP_PFO_X_i <- dplyr::recode(ISMP_PFO_items$ISMP_PFO_state, "BASE" = -0.5, "PAPE" = +0.5)


# simulate a sample of subjects
# calculate random intercept / random slope covariance
ISMP_PFO_covar <- ISMP_PFO_rho * ISMP_PFO_tau_0 * ISMP_PFO_tau_1

# put values into variance-covariance matrix
ISMP_PFO_cov_mx  <- matrix(
  c(ISMP_PFO_tau_0^2, ISMP_PFO_covar,
    ISMP_PFO_covar,   ISMP_PFO_tau_1^2),
  nrow = 2, byrow = TRUE)

# generate the by-subject random effects
ISMP_PFO_subject_rfx <- MASS::mvrnorm(n = ISMP_PFO_n_subj,
                                 mu = c(ISMP_PFO_T_0s = 0, ISMP_PFO_T_1s = 0),
                                 Sigma = ISMP_PFO_cov_mx)

# combine with subject IDs
ISMP_PFO_subjects <- data.frame(ISMP_PFO_subj_id = seq_len(ISMP_PFO_n_subj),
                                ISMP_PFO_subject_rfx)

# simulate a sample of subjects
# sample from a multivariate random distribution 
ISMP_PFO_subjects <- faux::rnorm_multi(
  n = ISMP_PFO_n_subj, 
  mu = 0, # means for random effects are always 0
  sd = c(ISMP_PFO_tau_0, ISMP_PFO_tau_1), # set SDs
  r = ISMP_PFO_rho, # set correlation, see ?faux::rnorm_multi
  varnames = c("ISMP_PFO_T_0s", "ISMP_PFO_T_1s")
)

# add subject IDs
ISMP_PFO_subjects$ISMP_PFO_subj_id <- seq_len(ISMP_PFO_n_subj)

# cross subject and item IDs; add an error term
# nrow(.) is the number of rows in the table
ISMP_PFO_trials <- crossing(ISMP_PFO_subjects, ISMP_PFO_items)  %>%
  mutate(ISMP_PFO_e_si = rnorm(nrow(.), mean = 0, sd = ISMP_PFO_sigma)) %>%
  dplyr::select(ISMP_PFO_subj_id, ISMP_PFO_item_id, ISMP_PFO_state, ISMP_PFO_X_i, everything())

ISMP_PFO_trials <- ISMP_PFO_trials[order(ISMP_PFO_trials$ISMP_PFO_subj_id, decreasing = F),]


# calculate the response variable
dat_sim_ISMP_PFO <- ISMP_PFO_trials %>%
  mutate(ISMP_PFO = ISMP_PFO_beta_0 + ISMP_PFO_T_0s + ISMP_PFO_O_0i + (ISMP_PFO_beta_1 + ISMP_PFO_T_1s) * ISMP_PFO_X_i + ISMP_PFO_e_si) %>%
  dplyr::select(ISMP_PFO_subj_id, ISMP_PFO_item_id, ISMP_PFO_state, ISMP_PFO_X_i, ISMP_PFO)

View(dat_sim_ISMP_PFO)



############### connect that into one simulated dataset ##########################

simulated_dataset_PAPE_isometric <- cbind(dat_sim_IBP_PFO[, c('IBP_PFO_subj_id',
                                               'IBP_PFO_item_id',
                                               'IBP_PFO_state',
                                               'IBP_PFO')], dat_sim_ISMP_PFO$ISMP_PFO)

colnames(simulated_dataset_PAPE_isometric)[colnames(simulated_dataset_PAPE_isometric) == "dat_sim_IBP_PFO$IBP_PFO"] ="IBP_PFO"
colnames(simulated_dataset_PAPE_isometric)[colnames(simulated_dataset_PAPE_isometric) == "dat_sim_ISMP_PFO$ISMP_PFO"] ="ISMP_PFO"
colnames(simulated_dataset_PAPE_isometric)[colnames(simulated_dataset_PAPE_isometric) == "IBP_PFO_state"] ="state"
colnames(simulated_dataset_PAPE_isometric)[colnames(simulated_dataset_PAPE_isometric) == "IBP_PFO_subj_id"] ="subj_ID"
colnames(simulated_dataset_PAPE_isometric)[colnames(simulated_dataset_PAPE_isometric) == "IBP_PFO_item_id"] ="item_ID"

simulated_dataset_PAPE_isometric$attempt_in_time_point <- rep(1:3)
simulated_dataset_PAPE_isometric$time_point <- rep(c(0, 5, 7, 10), each = 3)

unique_subj_ids <- unique(simulated_dataset_PAPE_isometric$subj_ID)

age_values <- rnorm(length(unique_subj_ids), mean = 28.42, sd = 7.0)
simulated_dataset_PAPE_isometric$age <- age_values[match(simulated_dataset_PAPE_isometric$subj_ID, unique_subj_ids)]

weight_values <- rnorm(length(unique_subj_ids), mean = 75.00, sd = 4.0)
simulated_dataset_PAPE_isometric$weight <- weight_values[match(simulated_dataset_PAPE_isometric$subj_ID, unique_subj_ids)]

height_values <- rnorm(length(unique_subj_ids), mean = 181.00, sd = 6.0)
simulated_dataset_PAPE_isometric$height <- height_values[match(simulated_dataset_PAPE_isometric$subj_ID, unique_subj_ids)]

IBP_one_rm_values <- rnorm(length(unique_subj_ids), mean = 90.00, sd = 10.0)
simulated_dataset_PAPE_isometric$IBP_one_rm <- IBP_one_rm_values[match(simulated_dataset_PAPE_isometric$subj_ID, unique_subj_ids)]
ISMP_one_rm_values <- rnorm(length(unique_subj_ids), mean = 50.00, sd = 5.0)
simulated_dataset_PAPE_isometric$ISMP_one_rm <- ISMP_one_rm_values[match(simulated_dataset_PAPE_isometric$subj_ID, unique_subj_ids)]

sex_values <- sample(c("M", "W"), length(unique_subj_ids), replace = TRUE)
simulated_dataset_PAPE_isometric$sex <- sex_values[match(simulated_dataset_PAPE_isometric$subj_ID, unique_subj_ids)]

simulated_dataset_PAPE_isometric <- simulated_dataset_PAPE_isometric[, c(1, 2, 13, 8:12, 3, 6, 7, 4:5)]

View(simulated_dataset_PAPE_isometric)


# -------------------------------------------- wrangling to wide format -------------------------------------- #

IBP_PFO_simulated_raw_data <- simulated_dataset_PAPE_isometric %>%
  unite("state_attempt", state, attempt_in_time_point, sep = "_") %>%
  pivot_wider(
    id_cols = c(subj_ID, sex, age, weight, height, IBP_one_rm, ISMP_one_rm),
    names_from = c("state_attempt", "time_point"),
    names_glue = "IBP_PFO_{state_attempt}_{time_point}",
    values_from = IBP_PFO 
  )

ISMP_PFO_simulated_raw_data <- simulated_dataset_PAPE_isometric %>%
  unite("state_attempt", state, attempt_in_time_point, sep = "_") %>%
  pivot_wider(
    id_cols = subj_ID,
    names_from = c("state_attempt", "time_point"),
    names_glue = "ISMP_PFO_{state_attempt}_{time_point}",
    values_from = ISMP_PFO 
  )

simulated_dataset_PAPE_isometric_wide <- cbind(IBP_PFO_simulated_raw_data, ISMP_PFO_simulated_raw_data[, c(2:25)])

setwd("E:/data/Statistics/Data")
#setwd("C:/Users/malir/OneDrive - Univerzita Karlova/Plocha/Data")

writexl::write_xlsx(simulated_dataset_PAPE_isometric, 'PAPE_isometric_simulated_dataset.xlsx')
writexl::write_xlsx(simulated_dataset_PAPE_isometric_wide, 'PAPE_isometric_simulated_dataset_wide.xlsx')
# writexl::write_xlsx(simulated_dataset_PAPE_isometric_wide_TP, 'PAPE_isometric_simulated_dataset_wide_TP.xlsx')

rm(list = ls())

simulated_dataset_PAPE_isometric <- readxl::read_excel("E:/data/Statistics/Data/PAPE_isometric_simulated_dataset.xlsx")
simulated_dataset_PAPE_isometric_wide <- readxl::read_excel("E:/data/Statistics/Data/PAPE_isometric_simulated_dataset_wide.xlsx")
# simulated_dataset_PAPE_isometric_wide_TP <- readxl::read_excel("E:/data/Statistics/Data/PAPE_isometric_simulated_dataset_wide_TP.xlsx")

# OR

simulated_dataset_PAPE_isometric <- readxl::read_excel("C:/Users/malir/OneDrive - Univerzita Karlova/Plocha/Data/PAPE_isometric_simulated_dataset.xlsx")
simulated_dataset_PAPE_isometric_wide <- readxl::read_excel("C:/Users/malir/OneDrive - Univerzita Karlova/Plocha/Data/PAPE_isometric_simulated_dataset_wide.xlsx")
# simulated_dataset_PAPE_isometric_wide_TP <- readxl::read_excel("C:/Users/malir/OneDrive - Univerzita Karlova/Plocha/Data/PAPE_isometric_simulated_dataset_wide_TP.xlsx")

# ------------------------------------------------------------------------------------------------------------ #
# ---------------------------------------------- descriptive stats ------------------------------------------- #
# ------------------------------------------------------------------------------------------------------------ #

descriptive <- round(pastecs::stat.desc(simulated_dataset_PAPE_isometric_wide[, c(3:55)], p = 0.95, norm = T), 2)
descriptive <- as.data.frame(descriptive)

(descriptive_body <- descriptive[c(4:6, 8:10, 13), c(1:53)])
print(descriptive_body)

(descriptive_body <- as.data.frame(t(descriptive_body)))
descriptive_body$variable <- rownames(descriptive_body)
descriptive_body <- descriptive_body[, c(8, 1:7)]
print(descriptive_body)

writexl::write_xlsx(descriptive_body, "descriptives_PAPE_isometric.xlsx")

# ------------------------------------------------------------------------------------------------------------ #
# ----------------------------------------- data visualization - rain cloud plots ---------------------------- #
# ------------------------------------------------------------------------------------------------------------ #

# ------------------------------------------ rain cloud plots based on test ---------------------------------- #

# simulated_dataset_PAPE_isometric$state <- as.factor(simulated_dataset_PAPE_isometric$state)

(simulated_dataset_PAPE_isometric_plot_IBP <- simulated_dataset_PAPE_isometric %>% 
    filter(state %in% c('BASE', 'PAPE')) %>% 
    ggplot(aes(x = factor(state), y = IBP_PFO, fill = factor(state))) +
    # add half-violin from {ggdist} package
    stat_halfeye(
      # adjust bandwidth
      adjust = 1.0,
      # move to the right
      justification = 0.0,
      # remove the slub interval
      .width = 0.0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = 0.12,
      # removing outliers
      outlier.color = NA,
      alpha = 0.5, position = position_dodge(width = 2.5)
    ) +
    stat_dots(
      # ploting on left side
      side = "left",
      # adjusting position
      justification = 1.1,
      # adjust grouping (binning) of observations
      binwidth = 0.35
    )+ theme(axis.line = element_line(colour = "black")) + theme_bw() + xlab("") + ylab("Peak force output (N)") + 
    scale_fill_manual(values = c("BASE" = "#AF46B4", "PAPE" = "#33a02c")))


(simulated_dataset_PAPE_isometric_plot_ISMP <- simulated_dataset_PAPE_isometric %>% 
    filter(state %in% c('BASE', 'PAPE')) %>% 
    ggplot(aes(x = factor(state), y = ISMP_PFO, fill = factor(state))) +
    # add half-violin from {ggdist} package
    stat_halfeye(
      # adjust bandwidth
      adjust = 1.0,
      # move to the right
      justification = 0.0,
      # remove the slub interval
      .width = 0.0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = 0.12,
      # removing outliers
      outlier.color = NA,
      alpha = 0.5, position = position_dodge(width = 2.5)
    ) +
    stat_dots(
      # ploting on left side
      side = "left",
      # adjusting position
      justification = 1.1,
      # adjust grouping (binning) of observations
      binwidth = 0.35
    )+ theme(axis.line = element_line(colour = "black")) + theme_bw() + xlab("") + ylab("Peak force output (N)") + 
    scale_fill_manual(values = c("BASE" = "#AF46B4", "PAPE" = "#33a02c")))


grid.arrange(simulated_dataset_PAPE_isometric_plot_IBP,
             simulated_dataset_PAPE_isometric_plot_ISMP, ncol=2)

# ------------------------------------ rain cloud plots based on time points -------------------------------- #

(simulated_dataset_PAPE_isometric_plot_TP_IBP_0 <- simulated_dataset_PAPE_isometric %>% 
    filter(state %in% c('BASE', 'PAPE'), time_point == 2) %>% 
    ggplot(aes(x = factor(state), y = IBP_PFO, fill = factor(state))) +
    # add half-violin from {ggdist} package
    stat_halfeye(
      # adjust bandwidth
      adjust = 1.0,
      # move to the right
      justification = 0.0,
      # remove the slub interval
      .width = 0.0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = 0.12,
      # removing outliers
      outlier.color = NA,
      alpha = 0.5, position = position_dodge(width = 2.5)
    ) +
    stat_dots(
      # ploting on left side
      side = "left",
      # adjusting position
      justification = 1.1,
      # adjust grouping (binning) of observations
      binwidth = 0.35
    )+ theme(axis.line = element_line(colour = "black")) + theme_bw() + xlab("IBP_PFO_2") + ylab("Peak force output (N)") + 
    scale_fill_manual(values = c("BASE" = "#AF46B4", "PAPE" = "#33a02c")))


(simulated_dataset_PAPE_isometric_plot_TP_IBP_5 <- simulated_dataset_PAPE_isometric %>% 
    filter(state %in% c('BASE', 'PAPE'), time_point == 5) %>% 
    ggplot(aes(x = factor(state), y = IBP_PFO, fill = factor(state))) +
    # add half-violin from {ggdist} package
    stat_halfeye(
      # adjust bandwidth
      adjust = 1.0,
      # move to the right
      justification = 0.0,
      # remove the slub interval
      .width = 0.0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = 0.12,
      # removing outliers
      outlier.color = NA,
      alpha = 0.5, position = position_dodge(width = 2.5)
    ) +
    stat_dots(
      # ploting on left side
      side = "left",
      # adjusting position
      justification = 1.1,
      # adjust grouping (binning) of observations
      binwidth = 0.35
    )+ theme(axis.line = element_line(colour = "black")) + theme_bw() + xlab("IBP_PFO_5") + ylab("Peak force output (N)") + 
    scale_fill_manual(values = c("BASE" = "#AF46B4", "PAPE" = "#33a02c")))


(simulated_dataset_PAPE_isometric_plot_TP_IBP_7 <- simulated_dataset_PAPE_isometric %>% 
    filter(state %in% c('BASE', 'PAPE'), time_point == 7) %>% 
    ggplot(aes(x = factor(state), y = IBP_PFO, fill = factor(state))) +
    # add half-violin from {ggdist} package
    stat_halfeye(
      # adjust bandwidth
      adjust = 1.0,
      # move to the right
      justification = 0.0,
      # remove the slub interval
      .width = 0.0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = 0.12,
      # removing outliers
      outlier.color = NA,
      alpha = 0.5, position = position_dodge(width = 2.5)
    ) +
    stat_dots(
      # ploting on left side
      side = "left",
      # adjusting position
      justification = 1.1,
      # adjust grouping (binning) of observations
      binwidth = 0.35
    )+ theme(axis.line = element_line(colour = "black")) + theme_bw() + xlab("IBP_PFO_7") + ylab("Peak force output (N)") + 
    scale_fill_manual(values = c("BASE" = "#AF46B4", "PAPE" = "#33a02c")))


(simulated_dataset_PAPE_isometric_plot_TP_IBP_10 <- simulated_dataset_PAPE_isometric %>% 
    filter(state %in% c('BASE', 'PAPE'), time_point == 10) %>% 
    ggplot(aes(x = factor(state), y = IBP_PFO, fill = factor(state))) +
    # add half-violin from {ggdist} package
    stat_halfeye(
      # adjust bandwidth
      adjust = 1.0,
      # move to the right
      justification = 0.0,
      # remove the slub interval
      .width = 0.0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = 0.12,
      # removing outliers
      outlier.color = NA,
      alpha = 0.5, position = position_dodge(width = 2.5)
    ) +
    stat_dots(
      # ploting on left side
      side = "left",
      # adjusting position
      justification = 1.1,
      # adjust grouping (binning) of observations
      binwidth = 0.35
    )+ theme(axis.line = element_line(colour = "black")) + theme_bw() + xlab("IBP_PFO_10") + ylab("Peak force output (N)") + 
    scale_fill_manual(values = c("BASE" = "#AF46B4", "PAPE" = "#33a02c")))


(simulated_dataset_PAPE_isometric_plot_TP_ISMP_0 <- simulated_dataset_PAPE_isometric %>% 
    filter(state %in% c('BASE', 'PAPE'), time_point == 2) %>% 
    ggplot(aes(x = factor(state), y = ISMP_PFO, fill = factor(state))) +
    # add half-violin from {ggdist} package
    stat_halfeye(
      # adjust bandwidth
      adjust = 1.0,
      # move to the right
      justification = 0.0,
      # remove the slub interval
      .width = 0.0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = 0.12,
      # removing outliers
      outlier.color = NA,
      alpha = 0.5, position = position_dodge(width = 2.5)
    ) +
    stat_dots(
      # ploting on left side
      side = "left",
      # adjusting position
      justification = 1.1,
      # adjust grouping (binning) of observations
      binwidth = 0.35
    )+ theme(axis.line = element_line(colour = "black")) + theme_bw() + xlab("ISMP_PFO_2") + ylab("Peak force output (N)") + 
    scale_fill_manual(values = c("BASE" = "#AF46B4", "PAPE" = "#33a02c")))


(simulated_dataset_PAPE_isometric_plot_TP_ISMP_5 <- simulated_dataset_PAPE_isometric %>% 
    filter(state %in% c('BASE', 'PAPE'), time_point == 5) %>% 
    ggplot(aes(x = factor(state), y = ISMP_PFO, fill = factor(state))) +
    # add half-violin from {ggdist} package
    stat_halfeye(
      # adjust bandwidth
      adjust = 1.0,
      # move to the right
      justification = 0.0,
      # remove the slub interval
      .width = 0.0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = 0.12,
      # removing outliers
      outlier.color = NA,
      alpha = 0.5, position = position_dodge(width = 2.5)
    ) +
    stat_dots(
      # ploting on left side
      side = "left",
      # adjusting position
      justification = 1.1,
      # adjust grouping (binning) of observations
      binwidth = 0.35
    )+ theme(axis.line = element_line(colour = "black")) + theme_bw() + xlab("ISMP_PFO_5") + ylab("Peak force output (N)") + 
    scale_fill_manual(values = c("BASE" = "#AF46B4", "PAPE" = "#33a02c")))


(simulated_dataset_PAPE_isometric_plot_TP_ISMP_7 <- simulated_dataset_PAPE_isometric %>% 
    filter(state %in% c('BASE', 'PAPE'), time_point == 7) %>% 
    ggplot(aes(x = factor(state), y = ISMP_PFO, fill = factor(state))) +
    # add half-violin from {ggdist} package
    stat_halfeye(
      # adjust bandwidth
      adjust = 1.0,
      # move to the right
      justification = 0.0,
      # remove the slub interval
      .width = 0.0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = 0.12,
      # removing outliers
      outlier.color = NA,
      alpha = 0.5, position = position_dodge(width = 2.5)
    ) +
    stat_dots(
      # ploting on left side
      side = "left",
      # adjusting position
      justification = 1.1,
      # adjust grouping (binning) of observations
      binwidth = 0.35
    )+ theme(axis.line = element_line(colour = "black")) + theme_bw() + xlab("ISMP_PFO_7") + ylab("Peak force output (N)") + 
    scale_fill_manual(values = c("BASE" = "#AF46B4", "PAPE" = "#33a02c")))


(simulated_dataset_PAPE_isometric_plot_TP_ISMP_10 <- simulated_dataset_PAPE_isometric %>% 
    filter(state %in% c('BASE', 'PAPE'), time_point == 10) %>% 
    ggplot(aes(x = factor(state), y = ISMP_PFO, fill = factor(state))) +
    # add half-violin from {ggdist} package
    stat_halfeye(
      # adjust bandwidth
      adjust = 1.0,
      # move to the right
      justification = 0.0,
      # remove the slub interval
      .width = 0.0,
      point_colour = NA
    ) +
    geom_boxplot(
      width = 0.12,
      # removing outliers
      outlier.color = NA,
      alpha = 0.5, position = position_dodge(width = 2.5)
    ) +
    stat_dots(
      # ploting on left side
      side = "left",
      # adjusting position
      justification = 1.1,
      # adjust grouping (binning) of observations
      binwidth = 0.35
    )+ theme(axis.line = element_line(colour = "black")) + theme_bw() + xlab("ISMP_PFO_10") + ylab("Peak force output (N)") + 
    scale_fill_manual(values = c("BASE" = "#AF46B4", "PAPE" = "#33a02c")))


grid.arrange(simulated_dataset_PAPE_isometric_plot_TP_IBP_0,
             simulated_dataset_PAPE_isometric_plot_TP_IBP_5,
             simulated_dataset_PAPE_isometric_plot_TP_IBP_7,
             simulated_dataset_PAPE_isometric_plot_TP_IBP_10,
             simulated_dataset_PAPE_isometric_plot_TP_ISMP_0,
             simulated_dataset_PAPE_isometric_plot_TP_ISMP_5,
             simulated_dataset_PAPE_isometric_plot_TP_ISMP_7,
             simulated_dataset_PAPE_isometric_plot_TP_ISMP_10, ncol=2)

# ------------------------------------------------------------------------------------------------------------ #
# ------------------------------------------------------ ICC ------------------------------------------------- #
# ------------------------------------------------------------------------------------------------------------ #

# PFO BASE and PAPE among all time points - IBP
cor_PFO_all_IBP <- corrplot(corr = cor(simulated_dataset_PAPE_isometric_wide[3:31], method = "spearman"), method="color",
                        bg = "grey10",
                        addgrid.col = "gray50", 
                        tl.cex=0.45,
                        order="hclust", 
                        number.cex=0.45,
                        addCoef.col = "black",
                        tl.col="black",
                        tl.srt=45,
                        sig.level = 0.05,
                        insig = "blank",
                        diag=TRUE,
                        col = colorRampPalette(c("yellow","green","brown1"))(100))

# PFO BASE and PAPE among all time points - ISMP
cor_PFO_all_ISMP <- corrplot(corr = cor(simulated_dataset_PAPE_isometric_wide[32:55], method = "spearman"), method="color",
                        bg = "grey10",
                        addgrid.col = "gray50", 
                        tl.cex=0.45,
                        order="hclust", 
                        number.cex=0.45,
                        addCoef.col = "black",
                        tl.col="black",
                        tl.srt=45,
                        sig.level = 0.05,
                        insig = "blank",
                        diag=TRUE,
                        col = colorRampPalette(c("yellow","green","brown1"))(100))

# ------------------ ICC between BASE and PAPE in different time points - IBP ------------------------- #

# PFO - ICC among attempts in time points - BASE
ICC_IBP_BASE_0 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('IBP_PFO_BASE_1_0', 'IBP_PFO_BASE_2_0', 'IBP_PFO_BASE_3_0')],
         model = "twoway", type = "consistency", unit = "average")
ICC_IBP_BASE_0_summary <- cbind("ICC value" = round(ICC_IBP_BASE_0$value, 3),
                               "Lower 95%CI" = round(ICC_IBP_BASE_0$lbound, 3),
                               "Upper 95%CI" = round(ICC_IBP_BASE_0$ubound, 3),
                               "p-value" = ICC_IBP_BASE_0$p.value)
rownames(ICC_IBP_BASE_0_summary) <- "ICC_IBP_BASE_0_summary"

ICC_IBP_BASE_5 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('IBP_PFO_BASE_1_5', 'IBP_PFO_BASE_2_5', 'IBP_PFO_BASE_3_5')],
         model = "twoway", type = "consistency", unit = "average")
ICC_IBP_BASE_5_summary <- cbind("ICC value" = round(ICC_IBP_BASE_5$value, 3),
                               "Lower 95%CI" = round(ICC_IBP_BASE_5$lbound, 3),
                               "Upper 95%CI" = round(ICC_IBP_BASE_5$ubound, 3),
                               "p-value" = ICC_IBP_BASE_5$p.value)
rownames(ICC_IBP_BASE_5_summary) <- "ICC_IBP_BASE_5_summary"

ICC_IBP_BASE_7 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('IBP_PFO_BASE_1_7', 'IBP_PFO_BASE_2_7', 'IBP_PFO_BASE_3_7')],
         model = "twoway", type = "consistency", unit = "average")
ICC_IBP_BASE_7_summary <- cbind("ICC value" = round(ICC_IBP_BASE_7$value, 3),
                                 "Lower 95%CI" = round(ICC_IBP_BASE_7$lbound, 3),
                                 "Upper 95%CI" = round(ICC_IBP_BASE_7$ubound, 3),
                                 "p-value" = ICC_IBP_BASE_7$p.value)
rownames(ICC_IBP_BASE_7_summary) <- "ICC_IBP_BASE_7_summary"

ICC_IBP_BASE_10 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('IBP_PFO_BASE_1_10', 'IBP_PFO_BASE_2_10', 'IBP_PFO_BASE_3_10')],
         model = "twoway", type = "consistency", unit = "average")
ICC_IBP_BASE_10_summary <- cbind("ICC value" = round(ICC_IBP_BASE_10$value, 3),
                                 "Lower 95%CI" = round(ICC_IBP_BASE_10$lbound, 3),
                                 "Upper 95%CI" = round(ICC_IBP_BASE_10$ubound, 3),
                                 "p-value" = ICC_IBP_BASE_10$p.value)
rownames(ICC_IBP_BASE_10_summary) <- "ICC_IBP_BASE_10_summary"

# PFO - ICC among attempts in time points - PAPE
ICC_IBP_PAPE_0 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('IBP_PFO_PAPE_1_0', 'IBP_PFO_PAPE_2_0', 'IBP_PFO_PAPE_3_0')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_IBP_PAPE_0_summary <- cbind("ICC value" = round(ICC_IBP_PAPE_0$value, 3),
                                 "Lower 95%CI" = round(ICC_IBP_PAPE_0$lbound, 3),
                                 "Upper 95%CI" = round(ICC_IBP_PAPE_0$ubound, 3),
                                 "p-value" = ICC_IBP_PAPE_0$p.value)
rownames(ICC_IBP_PAPE_0_summary) <- "ICC_IBP_PAPE_0_summary"

ICC_IBP_PAPE_5 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('IBP_PFO_PAPE_1_5', 'IBP_PFO_PAPE_2_5', 'IBP_PFO_PAPE_3_5')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_IBP_PAPE_5_summary <- cbind("ICC value" = round(ICC_IBP_PAPE_5$value, 3),
                                 "Lower 95%CI" = round(ICC_IBP_PAPE_5$lbound, 3),
                                 "Upper 95%CI" = round(ICC_IBP_PAPE_5$ubound, 3),
                                 "p-value" = ICC_IBP_PAPE_5$p.value)
rownames(ICC_IBP_PAPE_5_summary) <- "ICC_IBP_PAPE_5_summary"

ICC_IBP_PAPE_7 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('IBP_PFO_PAPE_1_7', 'IBP_PFO_PAPE_2_7', 'IBP_PFO_PAPE_3_7')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_IBP_PAPE_7_summary <- cbind("ICC value" = round(ICC_IBP_PAPE_7$value, 3),
                                 "Lower 95%CI" = round(ICC_IBP_PAPE_7$lbound, 3),
                                 "Upper 95%CI" = round(ICC_IBP_PAPE_7$ubound, 3),
                                 "p-value" = ICC_IBP_PAPE_7$p.value)
rownames(ICC_IBP_PAPE_7_summary) <- "ICC_IBP_PAPE_7_summary"

ICC_IBP_PAPE_10 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('IBP_PFO_PAPE_1_10', 'IBP_PFO_PAPE_2_10', 'IBP_PFO_PAPE_3_10')],
                            model = "twoway", type = "consistency", unit = "average")
ICC_IBP_PAPE_10_summary <- cbind("ICC value" = round(ICC_IBP_PAPE_10$value, 3),
                                  "Lower 95%CI" = round(ICC_IBP_PAPE_10$lbound, 3),
                                  "Upper 95%CI" = round(ICC_IBP_PAPE_10$ubound, 3),
                                  "p-value" = ICC_IBP_PAPE_10$p.value)
rownames(ICC_IBP_PAPE_10_summary) <- "ICC_IBP_PAPE_10_summary"

# -------------------- ICC between BASE and PAPE in different time points - ISMP ----------------------- #
# PFO - ICC among attempts in time points - BASE
ICC_ISMP_BASE_0 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('ISMP_PFO_BASE_1_0', 'ISMP_PFO_BASE_2_0', 'ISMP_PFO_BASE_3_0')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_ISMP_BASE_0_summary <- cbind("ICC value" = round(ICC_ISMP_BASE_0$value, 3),
                                "Lower 95%CI" = round(ICC_ISMP_BASE_0$lbound, 3),
                                "Upper 95%CI" = round(ICC_ISMP_BASE_0$ubound, 3),
                                "p-value" = ICC_ISMP_BASE_0$p.value)
rownames(ICC_ISMP_BASE_0_summary) <- "ICC_ISMP_BASE_0_summary"

ICC_ISMP_BASE_5 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('ISMP_PFO_BASE_1_5', 'ISMP_PFO_BASE_2_5', 'ISMP_PFO_BASE_3_5')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_ISMP_BASE_5_summary <- cbind("ICC value" = round(ICC_ISMP_BASE_5$value, 3),
                                "Lower 95%CI" = round(ICC_ISMP_BASE_5$lbound, 3),
                                "Upper 95%CI" = round(ICC_ISMP_BASE_5$ubound, 3),
                                "p-value" = ICC_ISMP_BASE_5$p.value)
rownames(ICC_ISMP_BASE_5_summary) <- "ICC_ISMP_BASE_5_summary"

ICC_ISMP_BASE_7 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('ISMP_PFO_BASE_1_7', 'ISMP_PFO_BASE_2_7', 'ISMP_PFO_BASE_3_7')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_ISMP_BASE_7_summary <- cbind("ICC value" = round(ICC_ISMP_BASE_7$value, 3),
                                "Lower 95%CI" = round(ICC_ISMP_BASE_7$lbound, 3),
                                "Upper 95%CI" = round(ICC_ISMP_BASE_7$ubound, 3),
                                "p-value" = ICC_ISMP_BASE_7$p.value)
rownames(ICC_ISMP_BASE_7_summary) <- "ICC_ISMP_BASE_7_summary"

ICC_ISMP_BASE_10 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('ISMP_PFO_BASE_1_10', 'ISMP_PFO_BASE_2_10', 'ISMP_PFO_BASE_3_10')],
                            model = "twoway", type = "consistency", unit = "average")
ICC_ISMP_BASE_10_summary <- cbind("ICC value" = round(ICC_ISMP_BASE_10$value, 3),
                                 "Lower 95%CI" = round(ICC_ISMP_BASE_10$lbound, 3),
                                 "Upper 95%CI" = round(ICC_ISMP_BASE_10$ubound, 3),
                                 "p-value" = ICC_ISMP_BASE_10$p.value)
rownames(ICC_ISMP_BASE_10_summary) <- "ICC_ISMP_BASE_10_summary"

# PFO - ICC among attempts in time points - PAPE
ICC_ISMP_PAPE_0 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('ISMP_PFO_PAPE_1_0', 'ISMP_PFO_PAPE_2_0', 'ISMP_PFO_PAPE_3_0')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_ISMP_PAPE_0_summary <- cbind("ICC value" = round(ICC_ISMP_PAPE_0$value, 3),
                                "Lower 95%CI" = round(ICC_ISMP_PAPE_0$lbound, 3),
                                "Upper 95%CI" = round(ICC_ISMP_PAPE_0$ubound, 3),
                                "p-value" = ICC_ISMP_PAPE_0$p.value)
rownames(ICC_ISMP_PAPE_0_summary) <- "ICC_ISMP_PAPE_0_summary"

ICC_ISMP_PAPE_5 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('ISMP_PFO_PAPE_1_5', 'ISMP_PFO_PAPE_2_5', 'ISMP_PFO_PAPE_3_5')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_ISMP_PAPE_5_summary <- cbind("ICC value" = round(ICC_ISMP_PAPE_5$value, 3),
                                "Lower 95%CI" = round(ICC_ISMP_PAPE_5$lbound, 3),
                                "Upper 95%CI" = round(ICC_ISMP_PAPE_5$ubound, 3),
                                "p-value" = ICC_ISMP_PAPE_5$p.value)
rownames(ICC_ISMP_PAPE_5_summary) <- "ICC_ISMP_PAPE_5_summary"

ICC_ISMP_PAPE_7 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('ISMP_PFO_PAPE_1_7', 'ISMP_PFO_PAPE_2_7', 'ISMP_PFO_PAPE_3_7')],
                           model = "twoway", type = "consistency", unit = "average")
ICC_ISMP_PAPE_7_summary <- cbind("ICC value" = round(ICC_ISMP_PAPE_7$value, 3),
                                "Lower 95%CI" = round(ICC_ISMP_PAPE_7$lbound, 3),
                                "Upper 95%CI" = round(ICC_ISMP_PAPE_7$ubound, 3),
                                "p-value" = ICC_ISMP_PAPE_7$p.value)
rownames(ICC_ISMP_PAPE_7_summary) <- "ICC_ISMP_PAPE_7_summary"

ICC_ISMP_PAPE_10 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('ISMP_PFO_PAPE_1_10', 'ISMP_PFO_PAPE_2_10', 'ISMP_PFO_PAPE_3_10')],
                            model = "twoway", type = "consistency", unit = "average")
ICC_ISMP_PAPE_10_summary <- cbind("ICC value" = round(ICC_ISMP_PAPE_10$value, 3),
                                 "Lower 95%CI" = round(ICC_ISMP_PAPE_10$lbound, 3),
                                 "Upper 95%CI" = round(ICC_ISMP_PAPE_10$ubound, 3),
                                 "p-value" = ICC_ISMP_PAPE_10$p.value)
rownames(ICC_ISMP_PAPE_10_summary) <- "ICC_ISMP_PAPE_10_summary"


# --------------------------------------- ICC - between BASE and PAPE ---------------------------------------- #

# ------------------- PFO - ICC between BASE and PAPE - IBP -------------------------- #
# PFO 2 BASE - PAPE
ICC_IBP_BASE_vs_PAPE_1_0 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('IBP_PFO_BASE_1_0', 'IBP_PFO_PAPE_1_0')],
         model = "twoway", type = "consistency", unit = "average")
ICC_IBP_BASE_vs_PAPE_1_0_summary <- cbind("ICC value" = round(ICC_IBP_BASE_vs_PAPE_1_0$value, 3),
                                  "Lower 95%CI" = round(ICC_IBP_BASE_vs_PAPE_1_0$lbound, 3),
                                  "Upper 95%CI" = round(ICC_IBP_BASE_vs_PAPE_1_0$ubound, 3),
                                  "p-value" = ICC_IBP_BASE_vs_PAPE_1_0$p.value)
rownames(ICC_IBP_BASE_vs_PAPE_1_0_summary) <- "ICC_IBP_BASE_vs_PAPE_1_0_summary"

ICC_IBP_BASE_vs_PAPE_2_0 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('IBP_PFO_BASE_2_0', 'IBP_PFO_PAPE_2_0')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_IBP_BASE_vs_PAPE_2_0_summary <- cbind("ICC value" = round(ICC_IBP_BASE_vs_PAPE_2_0$value, 3),
                                          "Lower 95%CI" = round(ICC_IBP_BASE_vs_PAPE_2_0$lbound, 3),
                                          "Upper 95%CI" = round(ICC_IBP_BASE_vs_PAPE_2_0$ubound, 3),
                                          "p-value" = ICC_IBP_BASE_vs_PAPE_2_0$p.value)
rownames(ICC_IBP_BASE_vs_PAPE_2_0_summary) <- "ICC_IBP_BASE_vs_PAPE_2_0_summary"

ICC_IBP_BASE_vs_PAPE_3_0 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('IBP_PFO_BASE_3_0', 'IBP_PFO_PAPE_3_0')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_IBP_BASE_vs_PAPE_3_0_summary <- cbind("ICC value" = round(ICC_IBP_BASE_vs_PAPE_3_0$value, 3),
                                          "Lower 95%CI" = round(ICC_IBP_BASE_vs_PAPE_3_0$lbound, 3),
                                          "Upper 95%CI" = round(ICC_IBP_BASE_vs_PAPE_3_0$ubound, 3),
                                          "p-value" = ICC_IBP_BASE_vs_PAPE_3_0$p.value)
rownames(ICC_IBP_BASE_vs_PAPE_3_0_summary) <- "ICC_IBP_BASE_vs_PAPE_3_0_summary"

# PFO 5 BASE - PAPE
ICC_IBP_BASE_vs_PAPE_1_5 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('IBP_PFO_BASE_1_5', 'IBP_PFO_PAPE_1_5')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_IBP_BASE_vs_PAPE_1_5_summary <- cbind("ICC value" = round(ICC_IBP_BASE_vs_PAPE_1_5$value, 3),
                                          "Lower 95%CI" = round(ICC_IBP_BASE_vs_PAPE_1_5$lbound, 3),
                                          "Upper 95%CI" = round(ICC_IBP_BASE_vs_PAPE_1_5$ubound, 3),
                                          "p-value" = ICC_IBP_BASE_vs_PAPE_1_5$p.value)
rownames(ICC_IBP_BASE_vs_PAPE_1_5_summary) <- "ICC_IBP_BASE_vs_PAPE_1_5_summary"

ICC_IBP_BASE_vs_PAPE_2_5 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('IBP_PFO_BASE_2_5', 'IBP_PFO_PAPE_2_5')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_IBP_BASE_vs_PAPE_2_5_summary <- cbind("ICC value" = round(ICC_IBP_BASE_vs_PAPE_2_5$value, 3),
                                          "Lower 95%CI" = round(ICC_IBP_BASE_vs_PAPE_2_5$lbound, 3),
                                          "Upper 95%CI" = round(ICC_IBP_BASE_vs_PAPE_2_5$ubound, 3),
                                          "p-value" = ICC_IBP_BASE_vs_PAPE_2_5$p.value)
rownames(ICC_IBP_BASE_vs_PAPE_2_5_summary) <- "ICC_IBP_BASE_vs_PAPE_2_5_summary"

ICC_IBP_BASE_vs_PAPE_3_5 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('IBP_PFO_BASE_3_5', 'IBP_PFO_PAPE_3_5')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_IBP_BASE_vs_PAPE_3_5_summary <- cbind("ICC value" = round(ICC_IBP_BASE_vs_PAPE_3_5$value, 3),
                                          "Lower 95%CI" = round(ICC_IBP_BASE_vs_PAPE_3_5$lbound, 3),
                                          "Upper 95%CI" = round(ICC_IBP_BASE_vs_PAPE_3_5$ubound, 3),
                                          "p-value" = ICC_IBP_BASE_vs_PAPE_3_5$p.value)
rownames(ICC_IBP_BASE_vs_PAPE_3_5_summary) <- "ICC_IBP_BASE_vs_PAPE_3_5_summary"

# PFO 7 BASE - PAPE
ICC_IBP_BASE_vs_PAPE_1_7 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('IBP_PFO_BASE_1_7', 'IBP_PFO_PAPE_1_7')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_IBP_BASE_vs_PAPE_1_7_summary <- cbind("ICC value" = round(ICC_IBP_BASE_vs_PAPE_1_7$value, 3),
                                          "Lower 95%CI" = round(ICC_IBP_BASE_vs_PAPE_1_7$lbound, 3),
                                          "Upper 95%CI" = round(ICC_IBP_BASE_vs_PAPE_1_7$ubound, 3),
                                          "p-value" = ICC_IBP_BASE_vs_PAPE_1_7$p.value)
rownames(ICC_IBP_BASE_vs_PAPE_1_7_summary) <- "ICC_IBP_BASE_vs_PAPE_1_7_summary"

ICC_IBP_BASE_vs_PAPE_2_7 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('IBP_PFO_BASE_2_7', 'IBP_PFO_PAPE_2_7')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_IBP_BASE_vs_PAPE_2_7_summary <- cbind("ICC value" = round(ICC_IBP_BASE_vs_PAPE_2_7$value, 3),
                                          "Lower 95%CI" = round(ICC_IBP_BASE_vs_PAPE_2_7$lbound, 3),
                                          "Upper 95%CI" = round(ICC_IBP_BASE_vs_PAPE_2_7$ubound, 3),
                                          "p-value" = ICC_IBP_BASE_vs_PAPE_2_2$p.value)
rownames(ICC_IBP_BASE_vs_PAPE_2_7_summary) <- "ICC_IBP_BASE_vs_PAPE_2_7_summary"

ICC_IBP_BASE_vs_PAPE_3_7 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('IBP_PFO_BASE_3_7', 'IBP_PFO_PAPE_3_7')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_IBP_BASE_vs_PAPE_3_7_summary <- cbind("ICC value" = round(ICC_IBP_BASE_vs_PAPE_3_7$value, 3),
                                          "Lower 95%CI" = round(ICC_IBP_BASE_vs_PAPE_3_7$lbound, 3),
                                          "Upper 95%CI" = round(ICC_IBP_BASE_vs_PAPE_3_7$ubound, 3),
                                          "p-value" = ICC_IBP_BASE_vs_PAPE_3_7$p.value)
rownames(ICC_IBP_BASE_vs_PAPE_3_7_summary) <- "ICC_IBP_BASE_vs_PAPE_3_7_summary"

# PFO 10 BASE - PAPE
ICC_IBP_BASE_vs_PAPE_1_10 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('IBP_PFO_BASE_1_10', 'IBP_PFO_PAPE_1_10')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_IBP_BASE_vs_PAPE_1_10_summary <- cbind("ICC value" = round(ICC_IBP_BASE_vs_PAPE_1_10$value, 3),
                                          "Lower 95%CI" = round(ICC_IBP_BASE_vs_PAPE_1_10$lbound, 3),
                                          "Upper 95%CI" = round(ICC_IBP_BASE_vs_PAPE_1_10$ubound, 3),
                                          "p-value" = ICC_IBP_BASE_vs_PAPE_1_10$p.value)
rownames(ICC_IBP_BASE_vs_PAPE_1_10_summary) <- "ICC_IBP_BASE_vs_PAPE_1_10_summary"

ICC_IBP_BASE_vs_PAPE_2_10 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('IBP_PFO_BASE_2_10', 'IBP_PFO_PAPE_2_10')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_IBP_BASE_vs_PAPE_2_10_summary <- cbind("ICC value" = round(ICC_IBP_BASE_vs_PAPE_2_10$value, 3),
                                          "Lower 95%CI" = round(ICC_IBP_BASE_vs_PAPE_2_10$lbound, 3),
                                          "Upper 95%CI" = round(ICC_IBP_BASE_vs_PAPE_2_10$ubound, 3),
                                          "p-value" = ICC_IBP_BASE_vs_PAPE_2_10$p.value)
rownames(ICC_IBP_BASE_vs_PAPE_2_10_summary) <- "ICC_IBP_BASE_vs_PAPE_2_10_summary"

ICC_IBP_BASE_vs_PAPE_3_10 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('IBP_PFO_BASE_3_10', 'IBP_PFO_PAPE_3_10')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_IBP_BASE_vs_PAPE_3_10_summary <- cbind("ICC value" = round(ICC_IBP_BASE_vs_PAPE_3_10$value, 3),
                                          "Lower 95%CI" = round(ICC_IBP_BASE_vs_PAPE_3_10$lbound, 3),
                                          "Upper 95%CI" = round(ICC_IBP_BASE_vs_PAPE_3_10$ubound, 3),
                                          "p-value" = ICC_IBP_BASE_vs_PAPE_3_10$p.value)
rownames(ICC_IBP_BASE_vs_PAPE_3_10_summary) <- "ICC_IBP_BASE_vs_PAPE_3_10_summary"

# ------------------- PFO - ICC between BASE and PAPE - ISMP ------------------------- #
# PFO 2 BASE - PAPE
ICC_ISMP_BASE_vs_PAPE_1_0 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('ISMP_PFO_BASE_1_0', 'ISMP_PFO_PAPE_1_0')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_ISMP_BASE_vs_PAPE_1_0_summary <- cbind("ICC value" = round(ICC_ISMP_BASE_vs_PAPE_1_0$value, 3),
                                          "Lower 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_1_0$lbound, 3),
                                          "Upper 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_1_0$ubound, 3),
                                          "p-value" = ICC_ISMP_BASE_vs_PAPE_1_0$p.value)
rownames(ICC_ISMP_BASE_vs_PAPE_1_0_summary) <- "ICC_ISMP_BASE_vs_PAPE_1_0_summary"

ICC_ISMP_BASE_vs_PAPE_2_0 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('ISMP_PFO_BASE_2_0', 'ISMP_PFO_PAPE_2_0')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_ISMP_BASE_vs_PAPE_2_0_summary <- cbind("ICC value" = round(ICC_ISMP_BASE_vs_PAPE_2_0$value, 3),
                                          "Lower 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_2_0$lbound, 3),
                                          "Upper 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_2_0$ubound, 3),
                                          "p-value" = ICC_ISMP_BASE_vs_PAPE_2_0$p.value)
rownames(ICC_ISMP_BASE_vs_PAPE_2_0_summary) <- "ICC_ISMP_BASE_vs_PAPE_2_0_summary"

ICC_ISMP_BASE_vs_PAPE_3_0 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('ISMP_PFO_BASE_3_0', 'ISMP_PFO_PAPE_3_0')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_ISMP_BASE_vs_PAPE_3_0_summary <- cbind("ICC value" = round(ICC_ISMP_BASE_vs_PAPE_3_0$value, 3),
                                          "Lower 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_3_0$lbound, 3),
                                          "Upper 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_3_0$ubound, 3),
                                          "p-value" = ICC_ISMP_BASE_vs_PAPE_3_0$p.value)
rownames(ICC_ISMP_BASE_vs_PAPE_3_0_summary) <- "ICC_ISMP_BASE_vs_PAPE_3_0_summary"

# PFO 5 BASE - PAPE
ICC_ISMP_BASE_vs_PAPE_1_5 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('ISMP_PFO_BASE_1_5', 'ISMP_PFO_PAPE_1_5')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_ISMP_BASE_vs_PAPE_1_5_summary <- cbind("ICC value" = round(ICC_ISMP_BASE_vs_PAPE_1_5$value, 3),
                                          "Lower 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_1_5$lbound, 3),
                                          "Upper 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_1_5$ubound, 3),
                                          "p-value" = ICC_ISMP_BASE_vs_PAPE_1_5$p.value)
rownames(ICC_ISMP_BASE_vs_PAPE_1_5_summary) <- "ICC_ISMP_BASE_vs_PAPE_1_5_summary"

ICC_ISMP_BASE_vs_PAPE_2_5 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('ISMP_PFO_BASE_2_5', 'ISMP_PFO_PAPE_2_5')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_ISMP_BASE_vs_PAPE_2_5_summary <- cbind("ICC value" = round(ICC_ISMP_BASE_vs_PAPE_2_5$value, 3),
                                          "Lower 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_2_5$lbound, 3),
                                          "Upper 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_2_5$ubound, 3),
                                          "p-value" = ICC_ISMP_BASE_vs_PAPE_2_5$p.value)
rownames(ICC_ISMP_BASE_vs_PAPE_2_5_summary) <- "ICC_ISMP_BASE_vs_PAPE_2_5_summary"

ICC_ISMP_BASE_vs_PAPE_3_5 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('ISMP_PFO_BASE_3_5', 'ISMP_PFO_PAPE_3_5')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_ISMP_BASE_vs_PAPE_3_5_summary <- cbind("ICC value" = round(ICC_ISMP_BASE_vs_PAPE_3_5$value, 3),
                                          "Lower 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_3_5$lbound, 3),
                                          "Upper 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_3_5$ubound, 3),
                                          "p-value" = ICC_ISMP_BASE_vs_PAPE_3_5$p.value)
rownames(ICC_ISMP_BASE_vs_PAPE_3_5_summary) <- "ICC_ISMP_BASE_vs_PAPE_3_5_summary"

# PFO 7 BASE - PAPE
ICC_ISMP_BASE_vs_PAPE_1_7 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('ISMP_PFO_BASE_1_7', 'ISMP_PFO_PAPE_1_7')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_ISMP_BASE_vs_PAPE_1_7_summary <- cbind("ICC value" = round(ICC_ISMP_BASE_vs_PAPE_1_7$value, 3),
                                          "Lower 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_1_7$lbound, 3),
                                          "Upper 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_1_7$ubound, 3),
                                          "p-value" = ICC_ISMP_BASE_vs_PAPE_1_7$p.value)
rownames(ICC_ISMP_BASE_vs_PAPE_1_7_summary) <- "ICC_ISMP_BASE_vs_PAPE_1_7_summary"

ICC_ISMP_BASE_vs_PAPE_2_7 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('ISMP_PFO_BASE_2_7', 'ISMP_PFO_PAPE_2_7')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_ISMP_BASE_vs_PAPE_2_7_summary <- cbind("ICC value" = round(ICC_ISMP_BASE_vs_PAPE_2_7$value, 3),
                                          "Lower 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_2_7$lbound, 3),
                                          "Upper 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_2_7$ubound, 3),
                                          "p-value" = ICC_ISMP_BASE_vs_PAPE_2_2$p.value)
rownames(ICC_ISMP_BASE_vs_PAPE_2_7_summary) <- "ICC_ISMP_BASE_vs_PAPE_2_7_summary"

ICC_ISMP_BASE_vs_PAPE_3_7 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('ISMP_PFO_BASE_3_7', 'ISMP_PFO_PAPE_3_7')],
                                     model = "twoway", type = "consistency", unit = "average")
ICC_ISMP_BASE_vs_PAPE_3_7_summary <- cbind("ICC value" = round(ICC_ISMP_BASE_vs_PAPE_3_7$value, 3),
                                          "Lower 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_3_7$lbound, 3),
                                          "Upper 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_3_7$ubound, 3),
                                          "p-value" = ICC_ISMP_BASE_vs_PAPE_3_7$p.value)
rownames(ICC_ISMP_BASE_vs_PAPE_3_7_summary) <- "ICC_ISMP_BASE_vs_PAPE_3_7_summary"

# PFO 10 BASE - PAPE
ICC_ISMP_BASE_vs_PAPE_1_10 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('ISMP_PFO_BASE_1_10', 'ISMP_PFO_PAPE_1_10')],
                                      model = "twoway", type = "consistency", unit = "average")
ICC_ISMP_BASE_vs_PAPE_1_10_summary <- cbind("ICC value" = round(ICC_ISMP_BASE_vs_PAPE_1_10$value, 3),
                                           "Lower 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_1_10$lbound, 3),
                                           "Upper 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_1_10$ubound, 3),
                                           "p-value" = ICC_ISMP_BASE_vs_PAPE_1_10$p.value)
rownames(ICC_ISMP_BASE_vs_PAPE_1_10_summary) <- "ICC_ISMP_BASE_vs_PAPE_1_10_summary"

ICC_ISMP_BASE_vs_PAPE_2_10 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('ISMP_PFO_BASE_2_10', 'ISMP_PFO_PAPE_2_10')],
                                      model = "twoway", type = "consistency", unit = "average")
ICC_ISMP_BASE_vs_PAPE_2_10_summary <- cbind("ICC value" = round(ICC_ISMP_BASE_vs_PAPE_2_10$value, 3),
                                           "Lower 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_2_10$lbound, 3),
                                           "Upper 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_2_10$ubound, 3),
                                           "p-value" = ICC_ISMP_BASE_vs_PAPE_2_10$p.value)
rownames(ICC_ISMP_BASE_vs_PAPE_2_10_summary) <- "ICC_ISMP_BASE_vs_PAPE_2_10_summary"

ICC_ISMP_BASE_vs_PAPE_3_10 <- irr::icc(simulated_dataset_PAPE_isometric_wide[, c('ISMP_PFO_BASE_3_10', 'ISMP_PFO_PAPE_3_10')],
                                      model = "twoway", type = "consistency", unit = "average")
ICC_ISMP_BASE_vs_PAPE_3_10_summary <- cbind("ICC value" = round(ICC_ISMP_BASE_vs_PAPE_3_10$value, 3),
                                           "Lower 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_3_10$lbound, 3),
                                           "Upper 95%CI" = round(ICC_ISMP_BASE_vs_PAPE_3_10$ubound, 3),
                                           "p-value" = ICC_ISMP_BASE_vs_PAPE_3_10$p.value)
rownames(ICC_ISMP_BASE_vs_PAPE_3_10_summary) <- "ICC_ISMP_BASE_vs_PAPE_3_10_summary"


ICC_all <- rbind(ICC_IBP_BASE_0_summary, ICC_IBP_BASE_5_summary, ICC_IBP_BASE_7_summary, ICC_IBP_BASE_10_summary,
      ICC_IBP_PAPE_0_summary, ICC_IBP_PAPE_5_summary, ICC_IBP_PAPE_7_summary, ICC_IBP_PAPE_10_summary,
      ICC_ISMP_BASE_0_summary, ICC_ISMP_BASE_5_summary, ICC_ISMP_BASE_7_summary, ICC_ISMP_BASE_10_summary,
      ICC_ISMP_PAPE_0_summary, ICC_ISMP_PAPE_5_summary, ICC_ISMP_PAPE_7_summary, ICC_ISMP_PAPE_10_summary,
      ICC_IBP_BASE_vs_PAPE_1_0_summary, ICC_IBP_BASE_vs_PAPE_2_0_summary, ICC_IBP_BASE_vs_PAPE_3_0_summary,
      ICC_IBP_BASE_vs_PAPE_1_5_summary, ICC_IBP_BASE_vs_PAPE_2_5_summary, ICC_IBP_BASE_vs_PAPE_3_5_summary,
      ICC_IBP_BASE_vs_PAPE_1_7_summary, ICC_IBP_BASE_vs_PAPE_2_7_summary, ICC_IBP_BASE_vs_PAPE_3_7_summary,
      ICC_IBP_BASE_vs_PAPE_1_10_summary, ICC_IBP_BASE_vs_PAPE_2_10_summary, ICC_IBP_BASE_vs_PAPE_3_10_summary,
      ICC_ISMP_BASE_vs_PAPE_1_0_summary, ICC_ISMP_BASE_vs_PAPE_2_0_summary, ICC_ISMP_BASE_vs_PAPE_3_0_summary,
      ICC_ISMP_BASE_vs_PAPE_1_5_summary, ICC_ISMP_BASE_vs_PAPE_2_5_summary, ICC_ISMP_BASE_vs_PAPE_3_5_summary,
      ICC_ISMP_BASE_vs_PAPE_1_7_summary, ICC_ISMP_BASE_vs_PAPE_2_7_summary, ICC_ISMP_BASE_vs_PAPE_3_7_summary,
      ICC_ISMP_BASE_vs_PAPE_1_10_summary, ICC_ISMP_BASE_vs_PAPE_2_10_summary, ICC_ISMP_BASE_vs_PAPE_3_10_summary)

ICC_all <- as.data.frame(ICC_all)

ICC_all$`p-value` <- ifelse(ICC_all$`p-value` <0.001, "<0.001", "<0.05") 
print(ICC_all)
range(ICC_all$`ICC value`)

# ----------------------------------------------------------------------------------------------------------- #
# ------------------------------------------------ mixed effect models -------------------------------------- #
# ----------------------------------------------------------------------------------------------------------- #

# ------------------------------------------------ IBP model ------------------------------------------------ #

# ------------------------ Response variable distribution -------------------- #

str(simulated_dataset_PAPE_isometric)

shapiro_test(simulated_dataset_PAPE_isometric$IBP_PFO)
ggplot(data = simulated_dataset_PAPE_isometric, aes(x = IBP_PFO)) + geom_density(fill = 'lightblue')

descdist(simulated_dataset_PAPE_isometric$IBP_PFO, boot = 1000)

# ------------------------------- Model settings ----------------------------- #

model_IBP_1 <- lme4::lmer(IBP_PFO ~ 1 + state + time_point + state:time_point + (1|subj_ID) + (1|subj_ID:attempt_in_time_point),
                          data = simulated_dataset_PAPE_isometric)

broom.mixed::tidy(model_IBP_1)

print(model_IBP_1)
summary(model_IBP_1)
(model_IBP_1_lmerTest <- lmerTest::lmer(IBP_PFO ~ 1 + state + time_point + state:time_point + (1|subj_ID) + (1|subj_ID:attempt_in_time_point),
                          data = simulated_dataset_PAPE_isometric))
summary(model_IBP_1_lmerTest)

# Variance estimation
VarCorr(model_IBP_1)

# R2M and R2C
# library(MuMIn)
r.squaredGLMM(model_IBP_1)

# Parameter estimates
# library(flexplot)
estimates(model_IBP_1)

# Visualize the model (sjPlot, flexplot packages)
# library(sjPlot)
sjPlot::plot_model(model_IBP_1, show.values = T, show.p = T)

# simulated_dataset_PAPE_isometric$state <- as.character(simulated_dataset_PAPE_isometric$state)
simulated_dataset_PAPE_isometric$state <- factor(simulated_dataset_PAPE_isometric$state, levels = c("BASE", "PAPE"), ordered = T)

(simulated_dataset_flexplot_IBP <- flexplot(IBP_PFO ~ state + time_point, 
                                            data = simulated_dataset_PAPE_isometric,
                                            model = model_IBP_1, jitter = NULL, spread = "stdev") + 
    theme(axis.text.x = element_text(angle = 0, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + 
    coord_cartesian(ylim = c(600, 1050)) + 
    scale_colour_manual(values = c("#4682B4", "#AF46B4", "#B47846", "#4BB446")) +
    scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) + 
    scale_shape_manual(values = c(19, 19, 19, 19)))

#################### and with 95%CIs based on the ggplot 2 ############################################################

summary_data_IBP <- simulated_dataset_PAPE_isometric %>%
  group_by(state, time_point) %>%
  summarise(mean_IBP = mean(IBP_PFO, na.rm = TRUE),
            sterr = sd(IBP_PFO, na.rm = TRUE) / sqrt(n()),
            lower_ci = mean_IBP - qt(0.975, df = n() - 1) * sterr,
            upper_ci = mean_IBP + qt(0.975, df = n() - 1) * sterr)

# Create the plot
(model_IBP_plot <- ggplot() +
    # 1. Individual participant data: gray points and connecting lines
    geom_line(data = simulated_dataset_PAPE_isometric, 
              aes(x = time_point, y = IBP_PFO, group = interaction(subj_ID, state), color = state), 
              alpha = 0.15) +  # Connect time points per subject, colored by state
    geom_point(data = simulated_dataset_PAPE_isometric, 
               aes(x = time_point, y = IBP_PFO, color = state), 
               size = 0.5, alpha = 0.2) +  # Individual data points
    
    # 2. Mean values with 95% confidence intervals
    geom_point(data = summary_data_IBP, 
               aes(x = time_point, y = mean_IBP, color = state), 
               size = 3, shape = 16) +  # Mean points
    geom_errorbar(data = summary_data_IBP, 
                  aes(x = time_point, ymin = lower_ci, ymax = upper_ci, color = state), 
                  width = 0.2, linewidth = 1.0) +  # 95% CI error bars
    geom_line(data = summary_data_IBP, 
              aes(x = time_point, y = mean_IBP, group = state, color = state), 
              linewidth = 1.0) +  # Line connecting means per state
    
    # Labels and theme
    xlab("Time Point") +
    ylab("Isometric Bench Press (N)") +
    scale_x_continuous(breaks = c(0, 5, 7, 10)) +  # Ensure time points appear correctly
    scale_color_manual(values = c("BASE" = "blue", "PAPE" = "red")) +  # Custom colors
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank()))  # Remove legend title

# ?pretty data frame? of predicted values for the response variable and its confidence interval
# library(ggeffects)
pred_model_IBP_1 <- ggpredict(model_IBP_1, terms = c("state"))
print(pred_model_IBP_1)

# --------------------------------- Multicollinearity -------------------------- #

# --------------- car package --------------- #

(VIF_model_IBP_1_car <- car::vif(model_IBP_1))

# ---------------------------- performance package ----------------------------- #

(VIF_model_IBP_1_performance <- performance::check_collinearity(model_IBP_1, component = "all"))
plot(VIF_model_IBP_1_performance)

VIF_model_IBP_1_performance <- as.data.frame(VIF_model_IBP_1_performance)
write_xlsx(VIF_model_IBP_1_performance, "E:/data/Statistics/Data/VIF_model_IBP_1_performance.xlsx")

# -------------------------------- VIF visualization --------------------------- #

par(mar = c(5, 10, 4, 8) + 0.2)

barplot(VIF_model_IBP_1_car, main = 'VIF values for the IBP model',beside = TRUE, horiz = T,
        col = 'orange', xlim = c(0,12), border = "black", axes = T,
        cex.names = 1.0, las = 1, xlab = "Variance Inflation Factor (car)")
abline(v = 10, col = "black", lty = 2)

# ----------------------------------- POST HOCS -------------------------------- #

# at first see if the interaction is significant
# (https://stats.stackexchange.com/questions/602206/posthoc-analysis-lmm-for-interaction-containing-a-factor-with-3-levels)

model_IBP_1_null <- lmerTest::lmer(IBP_PFO ~ state + time_point + (1|subj_ID) + (1|subj_ID:attempt_in_time_point),
                                   data = simulated_dataset_PAPE_isometric)

anova(model_IBP_1, model_IBP_1_null)

#library(emmeans)
lsmeans(model_IBP_1, pairwise ~ state:time_point)

simulated_dataset_PAPE_isometric$time_point <- as.character(simulated_dataset_PAPE_isometric$time_point) # for post hocks should be as character
model_IBP_1 <- lme4::lmer(IBP_PFO ~ 1 + state + time_point + state:time_point + (1|subj_ID) + (1|subj_ID:attempt_in_time_point),
                          data = simulated_dataset_PAPE_isometric)
str(simulated_dataset_PAPE_isometric)

(emm_IBP <- emmeans(model_IBP_1, specs = pairwise ~ state:time_point))

# ----------------------- 2 Homoscedasticity assumption ------------------------ #

#library(DHARMa)
res <- simulateResiduals(fittedModel = model_IBP_1, n = 10000)
residuals(res, quantileFunction = qnorm)

plot_res <- plotQQunif(res)
plot_res1 <- plotResiduals(res)

outliers(res, lowerQuantile = 0.025, upperQuantile = 0.975)

# plotting the residuals vs fitted values
plot(model_IBP_1, resid(., type = "pearson") ~ fitted(.), abline = 0,
     ylab = "Paerson's residuals", xlab = "Fitted values")
# there is no pattern within the residuals

# ----------------------- 4  Normal distribution of errors --------------------- #

# car package - qqPLot function
model_1_residuals <- resid(model_IBP_1)
qqPlot(model_1_residuals, pch = 20)

# ----------------------- 6 Outliers detection in the model -------------------- #

cooks_d_model_IBP_1 <- cooks.distance(model_IBP_1)
(cutoff_model_IBP_1 <- 4/(30))
plot(cooks_d_model_IBP_1, type = "b", ylab = "Cook's Distance",
     xlab = "Observation Index", pch = 20, col = "red", ylim = c(0, 0.15))
abline(h = cutoff_model_IBP_1,lty = 2)

# example with special cases for LMM/GLMM
# library(influence.ME)
estex_model_IBP_1 <- influence(model = model_IBP_1, "subj_ID")

dfbetas(estex_model_IBP_1)

plot(estex_model_IBP_1, which = "dfbetas", xlab="DFbetaS", ylab="Participant ID")

cooks.distance.estex(estex_model_IBP_1, sort = TRUE)

plot(estex_model_IBP_1, which = "cook",
     cutoff = cutoff_model_IBP_1, sort = TRUE,
     xlab = "Cook's Distance",
     ylab = "Participant ID", cex.axis = 0.5)

sigtest(estex = estex_model_IBP_1)


# ------------------------------------------------ ISMP model ----------------------------------------------- #

# ------------------------ Response variable distribution -------------------- #

str(simulated_dataset_PAPE_isometric)

shapiro_test(simulated_dataset_PAPE_isometric$ISMP_PFO)
ggplot(data = simulated_dataset_PAPE_isometric, aes(x = ISMP_PFO)) + geom_density(fill = 'lightblue')

descdist(simulated_dataset_PAPE_isometric$ISMP_PFO, boot = 1000)

# ------------------------------- Model settings ----------------------------- #

model_ISMP_1 <- lme4::lmer(ISMP_PFO ~ 1 + state + time_point + state:time_point + (1|subj_ID) + (1|subj_ID:attempt_in_time_point),
                          data = simulated_dataset_PAPE_isometric)

broom.mixed::tidy(model_ISMP_1)

print(model_ISMP_1)
summary(model_ISMP_1)
(model_ISMP_1_lmerTest <- lmerTest::lmer(ISMP_PFO ~ 1 + state + time_point + state:time_point + (1|subj_ID) + (1|subj_ID:attempt_in_time_point),
                                        data = simulated_dataset_PAPE_isometric))
summary(model_ISMP_1_lmerTest)

# Variance estimation
VarCorr(model_ISMP_1)

# R2M and R2C
# library(MuMIn)
r.squaredGLMM(model_ISMP_1)

# Parameter estimates
# library(flexplot)
estimates(model_ISMP_1)

# Visualize the model (sjPlot, flexplot packages)
# library(sjPlot)
sjPlot::plot_model(model_ISMP_1, show.values = T, show.p = T)

# simulated_dataset_PAPE_isometric$state <- as.character(simulated_dataset_PAPE_isometric$state)
simulated_dataset_PAPE_isometric$time_point <- as.numeric(simulated_dataset_PAPE_isometric$time_point) # go back to number again
simulated_dataset_PAPE_isometric$state <- factor(simulated_dataset_PAPE_isometric$state, levels = c("BASE", "PAPE"), ordered = T)

(simulated_dataset_flexplot_ISMP <- flexplot(ISMP_PFO ~ state + time_point, 
                                            data = simulated_dataset_PAPE_isometric,
                                            model = model_ISMP_1, jitter = NULL, spread = "stdev") + 
    theme(axis.text.x = element_text(angle = 0, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) + 
    coord_cartesian(ylim = c(100, 200)) + 
    scale_colour_manual(values = c("#4682B4", "#AF46B4", "#B47846", "#4BB446")) +
    scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) + 
    scale_shape_manual(values = c(19, 19, 19, 19)))


#################### and with 95%CIs based on the ggplot 2 ############################################################

summary_data_ISMP <- simulated_dataset_PAPE_isometric %>%
  group_by(state, time_point) %>%
  summarise(mean_ISMP = mean(ISMP_PFO, na.rm = TRUE),
            sterr = sd(ISMP_PFO, na.rm = TRUE) / sqrt(n()),
            lower_ci = mean_ISMP - qt(0.975, df = n() - 1) * sterr,
            upper_ci = mean_ISMP + qt(0.975, df = n() - 1) * sterr)
summary_data_IBP
# Create the plot
(model_ISMP_plot <- ggplot() +
    # 1. Individual participant data: gray points and connecting lines
    geom_line(data = simulated_dataset_PAPE_isometric, 
              aes(x = time_point, y = ISMP_PFO, group = interaction(subj_ID, state), color = state), 
              alpha = 0.15) +  # Connect time points per subject, colored by state
    geom_point(data = simulated_dataset_PAPE_isometric, 
               aes(x = time_point, y = ISMP_PFO, color = state), 
               size = 0.5, alpha = 0.2) +  # Individual data points
    
    # 2. Mean values with 95% confidence intervals
    geom_point(data = summary_data_ISMP, 
               aes(x = time_point, y = mean_ISMP, color = state), 
               size = 3, shape = 16) +  # Mean points
    geom_errorbar(data = summary_data_ISMP, 
                  aes(x = time_point, ymin = lower_ci, ymax = upper_ci, color = state), 
                  width = 0.2, linewidth = 1.0) +  # 95% CI error bars
    geom_line(data = summary_data_ISMP, 
              aes(x = time_point, y = mean_ISMP, group = state, color = state), 
              linewidth = 1.0) +  # Line connecting means per state
    
    # Labels and theme
    xlab("Time Point") +
    ylab("Isometric seated military press (N)") +
    scale_x_continuous(breaks = c(0, 5, 7, 10)) +  # Ensure time points appear correctly
    scale_color_manual(values = c("BASE" = "blue", "PAPE" = "red")) +  # Custom colors
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 0, hjust = 1), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank()))  # Remove legend title

# ?pretty data frame? of predicted values for the response variable and its confidence interval
# library(ggeffects)
pred_model_ISMP_1 <- ggpredict(model_ISMP_1, terms = c("state"))
print(pred_model_ISMP_1)

# --------------------------------- Multicollinearity -------------------------- #

# --------------- car package --------------- #

(VIF_model_ISMP_1_car <- car::vif(model_ISMP_1))

# ---------------------------- performance package ----------------------------- #

(VIF_model_ISMP_1_performance <- performance::check_collinearity(model_ISMP_1, component = "all"))
plot(VIF_model_ISMP_1_performance)

VIF_model_ISMP_1_performance <- as.data.frame(VIF_model_ISMP_1_performance)
write_xlsx(VIF_model_ISMP_1_performance, "E:/data/Statistics/Data/VIF_model_ISMP_1_performance.xlsx")

# -------------------------------- VIF visualization --------------------------- #

par(mar = c(5, 10, 4, 8) + 0.2)

barplot(VIF_model_ISMP_1_car, main = 'VIF values for the ISMP model',beside = TRUE, horiz = T,
        col = 'orange', xlim = c(0,12), border = "black", axes = T,
        cex.names = 1.0, las = 1, xlab = "Variance Inflation Factor (car)")
abline(v = 10, col = "black", lty = 2)

# ----------------------------------- POST HOCS -------------------------------- #

# at first see if the interaction is significant
# (https://stats.stackexchange.com/questions/602206/posthoc-analysis-lmm-for-interaction-containing-a-factor-with-3-levels)

model_ISMP_1_null <- lmerTest::lmer(ISMP_PFO ~ state + time_point + (1|subj_ID) + (1|subj_ID:attempt_in_time_point),
                                    data = simulated_dataset_PAPE_isometric)

anova(model_IBP_1, model_ISMP_1_null)

#library(emmeans)
lsmeans(model_ISMP_1, pairwise ~ state:time_point)

simulated_dataset_PAPE_isometric$time_point <- as.character(simulated_dataset_PAPE_isometric$time_point) # for post hocks should be as character
model_ISMP_1 <- lme4::lmer(ISMP_PFO ~ 1 + state + time_point + state:time_point + (1|subj_ID) + (1|subj_ID:attempt_in_time_point),
                          data = simulated_dataset_PAPE_isometric)
str(simulated_dataset_PAPE_isometric)

(emm <- emmeans(model_ISMP_1, specs = pairwise ~ state:time_point))

# ----------------------- 2 Homoscedasticity assumption ------------------------ #

#library(DHARMa)
res <- simulateResiduals(fittedModel = model_ISMP_1, n = 10000)
residuals(res, quantileFunction = qnorm)

plot_res <- plotQQunif(res)
plot_res1 <- plotResiduals(res)

outliers(res, lowerQuantile = 0.025, upperQuantile = 0.975)

# plotting the residuals vs fitted values
plot(model_ISMP_1, resid(., type = "pearson") ~ fitted(.), abline = 0,
     ylab = "Paerson's residuals", xlab = "Fitted values")
# there is no pattern within the residuals

# ----------------------- 4  Normal distribution of errors --------------------- #

# car package - qqPLot function
model_1_residuals <- resid(model_ISMP_1)
qqPlot(model_1_residuals, pch = 20)

# ----------------------- 6 Outliers detection in the model -------------------- #

cooks_d_model_ISMP_1 <- cooks.distance(model_ISMP_1)
(cutoff_model_ISMP_1 <- 4/(30))
plot(cooks_d_model_ISMP_1, type = "b", ylab = "Cook's Distance",
     xlab = "Observation Index", pch = 20, col = "red", ylim = c(0, 0.15))
abline(h = cutoff_model_ISMP_1,lty = 2)

# example with special cases for LMM/GLMM
# library(influence.ME)
estex_model_ISMP_1 <- influence(model = model_ISMP_1, "subj_ID")

dfbetas(estex_model_ISMP_1)

plot(estex_model_ISMP_1, which = "dfbetas", xlab="DFbetaS", ylab="Participant ID")

cooks.distance.estex(estex_model_ISMP_1, sort = TRUE)

plot(estex_model_ISMP_1, which = "cook",
     cutoff = cutoff_model_ISMP_1, sort = TRUE,
     xlab = "Cook's Distance",
     ylab = "Participant ID", cex.axis = 0.5)

sigtest(estex = estex_model_ISMP_1)

# ----------------------------------------------------------------------------------------------------------- #
# ------------------------------------------------ power analysis ------------------------------------------- #
# ----------------------------------------------------------------------------------------------------------- #

# ------------------- simulated data set function IBP ------------------- #

my_sim_data_IBP_PFO <- function(
    n_subj = 30,   # number of subjects
    n_BASE = 12,   # number of BASE
    n_PAPE = 12,   # number of PAPE
    beta_0 = 815,   # grand mean
    beta_1 = 76.5,   # effect of state
    omega_0 = 10,   # by-item random intercept sd
    tau_0 = 75,   # by-subject random intercept sd
    tau_1 = 30,   # by-subject random slope sd
    rho = 0.2,   # correlation between intercept and slope
    sigma = 20) { # residual (standard deviation)
  
  items <- data.frame(
    item_ID = seq_len(n_BASE + n_PAPE),
    state = rep(c("BASE", "PAPE"), c(n_BASE, n_PAPE)),
    X_i = rep(c(-0.5, 0.5), c(n_BASE, n_PAPE)),
    O_0i = rnorm(n = n_BASE + n_PAPE, mean = 0, sd = omega_0))
  
  # variance-covariance matrix
  cov_mx  <- matrix(
    c(tau_0^2, rho * tau_0 * tau_1, rho * tau_0 * tau_1, tau_1^2),
    nrow = 2, byrow = TRUE)
  
  subjects <- data.frame(subj_ID = seq_len(n_subj),
                         MASS::mvrnorm(n = n_subj,
                                       mu = c(T_0s = 0, T_1s = 0),
                                       Sigma = cov_mx))
  
  crossing(subjects, items) %>%
    mutate(e_si = rnorm(nrow(.), mean = 0, sd = sigma),
           IBP_PFO = beta_0 + T_0s + O_0i + (beta_1 + T_1s) * X_i + e_si,
           attempt_in_time_point = rep(1:3, length.out = nrow(.)),
           time_point = rep(c(0, 5, 7, 10), each = 3, length.out = nrow(.))) %>%
    dplyr::select(subj_ID, item_ID, state, X_i, IBP_PFO, attempt_in_time_point, time_point)
}

my_sim_data_IBP_PFO()

# ------------------- simulated data set function ISMP ------------------- #

my_sim_data_ISMP_PFO <- function(
    n_subj = 30,   # number of subjects
    n_BASE = 12,   # number of BASE
    n_PAPE = 12,   # number of PAPE
    beta_0 = 146,   # grand mean
    beta_1 = 10,   # effect of state
    omega_0 = 0,   # by-item random intercept sd
    tau_0 = 10,   # by-subject random intercept sd
    tau_1 = 4,   # by-subject random slope sd
    rho = 0.2,   # correlation between intercept and slope
    sigma = 3) { # residual (standard deviation)
  
  items <- data.frame(
    item_ID = seq_len(n_BASE + n_PAPE),
    state = rep(c("BASE", "PAPE"), c(n_BASE, n_PAPE)),
    X_i = rep(c(-0.5, 0.5), c(n_BASE, n_PAPE)),
    O_0i = rnorm(n = n_BASE + n_PAPE, mean = 0, sd = omega_0))
  
  # variance-covariance matrix
  cov_mx  <- matrix(
    c(tau_0^2, rho * tau_0 * tau_1, rho * tau_0 * tau_1, tau_1^2),
    nrow = 2, byrow = TRUE)
  
  subjects <- data.frame(subj_ID = seq_len(n_subj),
                         MASS::mvrnorm(n = n_subj,
                                       mu = c(T_0s = 0, T_1s = 0),
                                       Sigma = cov_mx))
  
  crossing(subjects, items) %>%
    mutate(e_si = rnorm(nrow(.), mean = 0, sd = sigma),
           ISMP_PFO = beta_0 + T_0s + O_0i + (beta_1 + T_1s) * X_i + e_si,
           attempt_in_time_point = rep(1:3, length.out = nrow(.)),
           time_point = rep(c(0, 5, 7, 10), each = 3, length.out = nrow(.))) %>%
    dplyr::select(subj_ID, item_ID, state, X_i, ISMP_PFO, attempt_in_time_point, time_point)
}

my_sim_data_ISMP_PFO()

# -------------------------------- single run function for power ---------------------------------- #

single_run_IBP_PFO <- function(...) {
  dat_sim <- my_sim_data_IBP_PFO(...)
  mod_sim <- lme4::lmer(IBP_PFO ~ 1 + state + time_point + state:time_point + (1|subj_ID) + (1|subj_ID:attempt_in_time_point),
                        data = dat_sim)
  table <- broom.mixed::tidy(mod_sim, effects = "fixed", conf.int = T)
  table <- as.data.frame(table)
  table <- table %>%
    dplyr::mutate(CI_includes_zero = ifelse(conf.low <= 0 & conf.high >= 0, 0, 1))
  table <- table %>%
    dplyr::mutate(CI_includes_zero_B = ifelse(conf.low <= 0 & conf.high >= 0, "Yes", "Not"))
  print(table)
}

single_run_IBP_PFO()

# --------------------------------------- number of simulations ---------------------------------- #

n_sims <- 10000

# ---------------------------------------- power IBP PFO ----------------------------------------- #

simulations_IBP_PFO <- purrr::map_df(1:n_sims, ~ single_run_IBP_PFO())

readr::write_csv(simulations_IBP_PFO, "E:/data/Statistics/Data/simulations_IBP_PFO.csv")

simulations_IBP_PFO <- readr::read_csv("E:/data/Statistics/Data/simulations_IBP_PFO.csv")

simulations_IBP_PFO %>% 
  filter(effect == "fixed") %>%
  group_by(term) %>%
  summarize(
    mean_estimate = mean(estimate),
    mean_se = mean(std.error),
    power = mean(CI_includes_zero > 0),
    power_B = mean(CI_includes_zero_B == "Not") * 100
  )

# ---------------------------------------- power simr IBP PFO ------------------------------------- #

str(simulated_dataset_PAPE_isometric)
simulated_dataset_PAPE_isometric$time_point <- as.numeric(simulated_dataset_PAPE_isometric$time_point)
model_IBP_1_simr <- lme4::lmer(IBP_PFO ~ 1 + state + time_point + state:time_point + (1|subj_ID) + (1|subj_ID:attempt_in_time_point),
                               data = simulated_dataset_PAPE_isometric)

power_model_1_1 <- powerSim(fit = model_IBP_1_simr, fixed('state'), nsim = 10)


print(power_model_1_1)
summary(power_model_1_1)

extended_model_1_1 <- simr::extend(model_IBP_1_simr, along = "subj_ID", n = 100)
pc_extended_model_1_1 <- powerCurve(extended_model_1_1,
                                    fixed('state'),
                                    along = 'subj_ID',
                                    nsim = 10, 
                                    alpha = 0.05, 
                                    progress = T)

print(summary(pc_extended_model_1_1))
plot(pc_extended_model_1_1)


# ----------------------------------------------------------------------------------------------------------- #
# ------------------------------------------------ equivalence testing -------------------------------------- #
# ---------------------------------------------------- TOSTER ----------------------------------------------- #

# ----------------------------------------------------------------------------------------------------------- #

# ----------------------------- single run function for equivalence ------------------------------- #

single_run_equivalance <- function(...){
  dat_IBP <- my_sim_data_IBP_PFO(...) %>%
    mutate(IBP_PFO = scale(IBP_PFO, center = TRUE, scale = TRUE))
  dat_ISMP <- my_sim_data_ISMP_PFO(...) %>%
    mutate(ISMP_PFO = scale(ISMP_PFO, center = TRUE, scale = TRUE))
  
  # Models fitting
  dat_IBP$time_point <- as.character(dat_IBP$time_point)
  dat_ISMP$time_point <- as.character(dat_ISMP$time_point)
  (mod_IBP <- lme4::lmer(IBP_PFO ~ 1 + state + time_point + state:time_point + (1|subj_ID) + (1|subj_ID:attempt_in_time_point),
                         data = dat_IBP))
  (mod_ISMP <- lme4::lmer(ISMP_PFO ~ 1 + state + time_point + state:time_point + (1|subj_ID) + (1|subj_ID:attempt_in_time_point),
                          data = dat_ISMP))
  
  # pairwise comparisons to extract the exact time points differences
  # (emm_IBP <- emmeans(mod_IBP, specs = pairwise ~ state:time_point))
  # (emm_ISMP <- emmeans(mod_ISMP, specs = pairwise ~ state:time_point))
  
  # Extracting the fixed effects
  fixef_IBP <- broom.mixed::tidy(mod_IBP, effects = "fixed") %>% filter(term == "statePAPE")
  fixef_ISMP <- broom.mixed::tidy(mod_ISMP, effects = "fixed") %>% filter(term == "statePAPE")
  
  # extracting effect sizes
  d_mod1 <- fixef_IBP$estimate #/ sigma_1)
  d_mod2 <- fixef_ISMP$estimate #/ sigma_2)
  #(d_mod1_sig <- fixef_IBP$estimate/ sigma_1)
  #(d_mod2_sig <- fixef_ISMP$estimate/ sigma_2)
  
  # Extract estimates and standard errors
  # (estimate_IBP <- fixef_IBP$estimate)
  # (estimate_ISMP <- fixef_ISMP$estimate)
  # (se_IBP <- fixef_IBP$std.error)
  # (se_ISMP <- fixef_ISMP$std.error)
  sigma_1 <- sigma(mod_IBP)
  sigma_2 <- sigma(mod_ISMP)
  
  low_eqbound <- -0.5
  high_eqbound <- 0.5
  
  tost_result <- TOSTER::tsum_TOST(
    n1 = 30, 
    n2 = 30,
    m1 = d_mod1,
    m2 = d_mod2,
    sd1 = sigma_1,
    sd2 = sigma_2,
    r12 = 0.75, # 0.75 better
    low_eqbound = low_eqbound,
    high_eqbound = high_eqbound,
    paired = T
  )
  table <- as.data.frame(tost_result$effsize)
  table <- table %>%
    dplyr::mutate(CI_includes_zero = ifelse(lower.ci <= 0 & upper.ci >= 0, 0, 1))
  table <- table %>%
    dplyr::mutate(CI_includes_zero_B = ifelse(lower.ci <= 0 & upper.ci >= 0, "Yes", "Not"))
  print(table)
}

single_run_equivalance()

# --------------------------------- power of equivalence --------------------------------------- #

simulations_PAPE_isometric_equivalence <- purrr::map_df(1:n_sims, ~ single_run_equivalance())

readr::write_csv(simulations_PAPE_isometric_equivalence, "E:/data/Statistics/Data/simulations_PAPE_isometric_equivalence.csv")

simulations_PAPE_isometric_equivalence <- readr::read_csv("E:/data/Statistics/Data/simulations_PAPE_isometric_equivalence.csv")

simulations_PAPE_isometric_equivalence %>% 
  summarize(
    mean_estimate = mean(estimate),
    mean_se = mean(SE),
    power = mean(CI_includes_zero > 0),
    power_B = mean(CI_includes_zero_B == "Not") * 100
  )

# ----------------------------------------------------------------------------------------------- #

