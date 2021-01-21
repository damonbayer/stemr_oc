library(rms)
library(tidyverse)
library(lubridate)

# Setup -------------------------------------------------------------------
cases_age_lag <- 21
deaths_csv <- "~/Documents/oc_covid19/data/oc/12.7.20 release to UCI team.csv"
neg_line_list_csv <- "~/Documents/oc_covid19/data/oc/All ELR PCR tests updated 12.7.20.csv"

source("~/Documents/oc_covid19/data/oc/oc_read_data.R")

# IFR Model ---------------------------------------------------------------
nature_dat <-  tibble(age = c(2, 7, 12, 17, 22, 27, 32, 37, 42, 47, 52, 57, 62, 67, 72, 77, 85),
                      ifr = c(3e-05, 1e-05, 1e-05, 3e-05, 6e-05, 0.00013, 0.00024, 4e-04,
                              0.00075, 0.00121, 0.00207, 0.00323, 0.00456, 0.01075, 0.01674,
                              0.03203, 0.08292))

nature_dat_model <- glm(formula = ifr ~ age, family = gaussian(link = "logit"), data = nature_dat)

tibble(age = seq(0,90,5)) %>%
  mutate(., ifr = predict(nature_dat_model, newdata = ., type = "response")) %>%
  ggplot(aes(age, ifr)) +
  geom_line() +
  geom_point() +
  scale_y_log10(label = scales::percent)


tibble()

351/2979

seroprev_data <-

# Process Data ------------------------------------------------------------
process_and_save_data <- function(snf_omit) {
  new_deaths_tbl <-
    read_csv(deaths_csv,
             col_types = cols(.default = col_skip(),
                              `HighRisk` = col_character(),
                              `DtDeath` = col_date("%Y-%m-%d"),
                              `DeathDueCOVID` = col_character())) %>%
    replace_na(list(HighRisk = "No")) %>%
    filter(DeathDueCOVID == "Y") %>%
    purrr::when(snf_omit ~ filter(., !str_detect(HighRisk, "Resident")),
                ~ .) %>%
    select(date = DtDeath) %>%
    count(date, name = "deaths")

  neg_line_list <- read_csv(neg_line_list_csv,
                            col_types = cols(.default = col_skip(),
                                             Age = col_double(),
                                             PersonId = col_integer(),
                                             Specimen.Collected.Date = col_date("%m-%d-%Y"),
                                             Resulted.Organism = col_character())) %>%
    filter(!is.na(Resulted.Organism)) %>%
    mutate(test_result = fct_collapse(str_to_lower(Resulted.Organism),
                                      negative = negative_test_synonyms,
                                      positive = positive_test_synonyms,
                                      other = other_test_synonyms)) %>%
    select(id = PersonId, age = Age, date = Specimen.Collected.Date, test_result) %>%
    filter(date >= lubridate::ymd("2020-01-01")) %>%
    group_by(id) %>%
    arrange(date) %>%
    ungroup()

  if(length(levels(neg_line_list$test_result)) != 3) stop("New test result category not accounted for.")

  first_pos <- neg_line_list %>%
    filter(test_result == "positive") %>%
    group_by(id) %>%
    summarise(first_pos = min(date))

  neg_line_list_filtered <- left_join(neg_line_list, first_pos) %>%
    mutate(first_pos = replace_na(first_pos, lubridate::ymd("9999-12-31"))) %>%
    filter(date <= first_pos) %>%
    select(-first_pos) %>%
    distinct()

  dat_for_ifr_model <- neg_line_list %>%
    filter(test_result == "positive") %>%
    filter(date >= "2020-03-14",
           date <= max(neg_line_list[["date"]]) - 14) %>%
    select(date, age) %>%
    mutate(., ifr = predict(nature_dat_model, newdata = ., type = "response")) %>%
    mutate(date = as.numeric(date))

  rcs_model_ifr <- ols(ifr ~ rcs(date,
                                 # quantile(date, c(.05, .23, .41, .59, .77, 0.95))),
                                 quantile(date, c(.025, .1833, .3417, .5, .6583, .8167, .975))),
                       data = dat_for_ifr_model %>% drop_na())


  final_dat <- neg_line_list_filtered %>%
    count(date, test_result) %>%
    pivot_wider(names_from = test_result, values_from = n) %>%
    full_join(new_deaths_tbl) %>%
    right_join(., tibble(date = seq(min(.[["date"]]), max(.[["date"]]), by = 1))) %>%
    arrange(date) %>%
    replace(is.na(.), 0) %>%
    mutate(cases = positive,
           tests = negative + positive + other,
           lagged_pos_ifr = unname(predict(rcs_model_ifr,
                                           newdata = tibble(date = as.numeric(date) - cases_age_lag)))) %>%
    select(date, deaths, cases, tests, lagged_pos_ifr)

  ggplot(final_dat, aes(date, lagged_pos_ifr)) +
    geom_line() +
    scale_x_date(limits = c(ymd("2020-03-01"), max(final_dat$date)))

  if(snf_omit) write_rds(final_dat, "data/oc_data_no_snf.rds")
  if(!snf_omit) write_rds(final_dat, "data/oc_data.rds")

}

process_and_save_data(snf_omit = T)
process_and_save_data(snf_omit = F)
