# generate the phase 1 dataset for the HVTN505 analysis
library("dplyr")
library("tidyr")
library("tibble")
library("readr")
# ------------------------------------------------------
# Read in the full dataset and make sure that
# strata are computed across ALL participants
# ------------------------------------------------------
full_data <- readr::read_csv("data/primary505_for_sharing.csv") 
full_data_with_strata <- full_data %>% 
  mutate(hg_strata = ifelse(is.na(hg_strata), 
                            paste0(
                              ifelse(trt == 0, "Plac/", "Vacc/"),
                              ifelse(racecc == "White", "White", "Blk_Hisp"),
                              ifelse(trt == 1 & is.na(BMI), "/Missing", ""),
                              ifelse(trt == 1 & BMI < 25 & !is.na(BMI), "/[18.4, 25)", ""),
                              ifelse(trt == 1 & 25 <= BMI & BMI < 29.8 & !is.na(BMI), "/[25, 29.8)", ""),
                              ifelse(trt == 1 & 29.8 <= BMI & !is.na(BMI), "/[29.8, 40)", "")
                            ), 
                            hg_strata)) %>% 
  filter(!is.na(BMI))

# read in table to match pub_id from this dataset with the
# ID in the HVTN505 R package
pub_id_converter <- readr::read_csv("data/rx_v2_bdw.csv") %>% 
  select(hvtn505_id, primary505_id) %>% 
  mutate(hvtn505_id = gsub("-", "", hvtn505_id))

# update pub id
full_data_with_matched_ids <- full_data_with_strata %>% 
  left_join(pub_id_converter %>% rename(pub_id = primary505_id), by = "pub_id")
# ------------------------------------------------------
# Compute inverse probability weights
# ------------------------------------------------------
# compute strata 
full_data_with_stratuminds <- full_data_with_matched_ids %>% 
  mutate(
    stratuminds = case_when(
      hg_strata == "Plac/Blk_Hisp" ~ 1,
      hg_strata == "Plac/White" ~ 2,
      hg_strata == "Vacc/Blk_Hisp/[18.4, 25)" ~ 3,
      hg_strata == "Vacc/Blk_Hisp/[25, 29.8)" ~ 4,
      hg_strata == "Vacc/Blk_Hisp/[29.8, 40)" ~ 5,
      hg_strata == "Vacc/White/[18.4, 25)" ~ 6,
      hg_strata == "Vacc/White/[25, 29.8)" ~ 7,
      hg_strata == "Vacc/White/[29.8, 40)" ~ 8
    ),
    stratuminds_vaccs = stratuminds - 2
  )
cc_data <- full_data_with_stratuminds %>% 
  filter(casecontrol == 1)
# numbers of cases and controls in each stratum
n0 <- full_data_with_stratuminds %>% 
  filter(cc_cohort == 1, HIVwk28preunbl == 0) %>% 
  group_by(cc_strata) %>% 
  summarize(n0 = n(), .groups = "drop")
n1 <- full_data_with_stratuminds %>% 
  filter(cc_cohort == 1, HIVwk28preunbl == 1) %>% 
  group_by(cc_strata) %>% 
  summarize(n1 = n(), .groups = "drop")
# compute the weights
weights <- cbind(n0$n0, n1$n1) / table(cc_data$stratuminds, cc_data$HIVwk28preunbl)
full_data_with_weights <- full_data_with_stratuminds %>% 
  mutate(weight = ifelse(HIVwk28preunbl == 1, 
                         1, 
                         weights[, 1][full_data_with_stratuminds$stratuminds]))

# -----------------------------------------------------
# Save off the tibble with weights and ID, Z covariates
# -----------------------------------------------------
Z_plus_weights <- full_data_with_weights %>% 
  select(hvtn505_id, trt, HIVpreunbl, age, BMI, bhvrisk, weight) %>% 
  rename(ptid = hvtn505_id, Y = HIVpreunbl)

# first check to make sure that things are correct
library("HVTN505")
data("dat.505")
sum(dat.505$ptid %in% Z_plus_weights$ptid)
# are the Z covariates the same? Answer: yes
Z_plus_weights %>% 
  filter(ptid %in% dat.505$ptid) %>% 
  select(-weight, -trt)
head(dat.505[, c("ptid", "case", "age", "BMI", "bhvrisk")], 10)
mean(Z_plus_weights %>% 
      filter(ptid %in% dat.505$ptid) %>% 
       select(-weight, -trt) == dat.505[, c("ptid", "case", "age", "BMI", "bhvrisk")])
# Are the weights the same? Answer: yes (to 13 digits, which is close enough)
Z_plus_weights %>% 
  filter(ptid %in% dat.505$ptid) %>% 
  select(ptid, weight)
head(dat.505[, c("ptid", "wt")], 10)
Z_plus_weights %>% 
  filter(ptid %in% dat.505$ptid) %>% 
  select(ptid, weight) %>% 
  mutate(ptid = as.numeric(ptid)) %>% 
  left_join(dat.505 %>% select(ptid, wt), by = "ptid") %>% 
  print(n = Inf)
mean(Z_plus_weights %>% 
       filter(ptid %in% dat.505$ptid) %>% 
       select(ptid, weight)  == dat.505[, c("ptid", "wt")])
dgts <- 13
mean(Z_plus_weights %>% 
       filter(ptid %in% dat.505$ptid) %>% 
       select(ptid, weight) %>% pull(weight) %>% round(digits = dgts) 
     == dat.505[, c("ptid", "wt")] %>% pull(wt) %>%  round(digits = dgts))

# minor cleanup prior to saving
Z_plus_weights_save <- Z_plus_weights %>% 
  mutate(ptid = as.numeric(ptid)) %>% 
  filter(!is.na(Y), !is.na(age), !is.na(BMI), !is.na(bhvrisk))

saveRDS(Z_plus_weights_save, file = "data/z_and_weights_for_505_analysis.rds")
