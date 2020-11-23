library(rms)
library(tidyverse)
library(lubridate)

# Setup -------------------------------------------------------------------
snf_omit <- F
cases_age_lag <- 18
deaths_csv <- "data/11.9.20 release to UCI team.csv"
neg_line_list_csv <- "data/All ELR PCR tests updated 11.9.20.csv"

negative_test_synonyms <- c("not detected",
                            "negative",
                            "coronavirus 2019 novel not detected",
                            "negative (qualifier value)",
                            "not detected (qualifier value)",
                            "sars cov-2 negative",
                            "undetected",
                            "inst_negative",
                            "neg-see report",
                            "sars-cov-2 rna not detected by naa",
                            "none detected",
                            "not detected in pooled specimen",
                            "not detected in pooled specimen (qualifier value)")

positive_test_synonyms <- c("detected",
                            "coronavirus 2019 novel positive",
                            "positive",
                            "positive (qualifier value)",
                            "sars cov-2 positive",
                            "detected (qualifier value)",
                            "presumptive pos",
                            "positive for 2019-ncov",
                            "presumptive positive",
                            "coronavirus 2019 novel presumptive pos",
                            "coronavirus 2019 novel detected",
                            "yes",
                            "coronavirus 2019 novel",
                            "presumptive positive for 2019-ncov",
                            "sars cov-2 presumptive pos",
                            "presumptive pos. for 2019-ncov",
                            "presumptive positive (qualifier value)",
                            "presumptive detected",
                            "reactive",
                            "sars-cov-2",
                            "interpretive information: 2019 novel coronavirus sars-cov-2 by pcr")

other_test_synonyms <- c("inconclusive",
                         "indeterminate",
                         "specimen unsatisfactory",
                         "invalid",
                         "test not performed",
                         "not provided (qualifier value)",
                         "see comment",
                         "tnp",
                         "coronavirus 2019 novel inconclusive",
                         "not tested",
                         "phoned results (and readback confirmed) to:",
                         "see note",
                         "clotted",
                         "coronavirus 2019 novel unsatisfactory",
                         # "cryptococcus neoformans",
                         "equivocal",
                         "non reactive",
                         "result comments",
                         "sars cov-2 inconclusive",
                         "test not done",
                         "test not perf",
                         "not pregnant",
                         "biofiresarsneg",
                         "equivocal result",
                         "coronavirus 2019 novel inconcluside",
                         "unsatisfactory",
                         "undefined",
                         "*corrected report* by",
                         "specimen unsatifactory for evaluation",
                         "warning....please disregard results.",
                         "presumptive result to be confirmed",
                         "indeterminate (qualifier value)",
                         "invalid result",
                         "specimen unsatisfactory for evaluation",
                         "specimen received mislabeled",
                         "enterococcus faecalis",
                         "carbapenem resistant pseudomonas aeruginosa",
                         "enterobacter cloacae complex (organism)",
                         "unknown",
                         "multiple drug-resistant serratia marcescens",
                         "genus enterococcus",
                         "acinetobacter baumannii (organism)"
)


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



# Read Data ---------------------------------------------------------------
new_deaths_tbl <-
  read_csv(deaths_csv,
           col_types = cols(.default = col_skip(),
                            `HighRisk` = col_character(),
                            `DtDeath` = col_date("%Y-%m-%d"),
                            `DeathDueCOVID` = col_character())) %>%
  replace_na(list(HighRisk = "No")) %>%
  filter(DeathDueCOVID == "Y") %>%
  when(snf_omit ~ filter(., !str_detect(HighRisk, "Resident")),
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


# Process Data ------------------------------------------------------------
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


ggplot(final_dat, aes(date, lagged_pos_ifr)) +
  geom_line() +
  scale_x_date(limits = c(ymd("2020-03-01"), max(final_dat$date)))


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

if(snf_omit) write_rds(final_dat, "data/oc_data_no_snf.rds")
if(!snf_omit) write_rds(final_dat, "data/oc_data.rds")

