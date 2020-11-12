library(rms)
library(tidyverse)
cases_age_lag <- 14

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

new_deaths_tbl <- read_csv("data/11.9.20 release to UCI team.csv",
         col_types = cols(.default = col_skip(),
                          `DtDeath` = col_date("%Y-%m-%d"),
                          `DeathDueCOVID` = col_character())) %>%
  filter(DeathDueCOVID == "Y") %>%
  select(date = DtDeath) %>%
  count(date, name = "deaths")


neg_line_list <- read_csv("data/All ELR PCR tests updated 11.9.20.csv",
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

if(length(levels(neg_line_list$test_result)) != 3) warning("New test result category not accounted for.")

first_pos <- neg_line_list %>%
  filter(test_result == "positive") %>%
  group_by(id) %>%
  summarise(first_pos = min(date))

neg_line_list_filtered <- left_join(neg_line_list, first_pos) %>%
  mutate(first_pos = replace_na(first_pos, lubridate::ymd("9999-12-31"))) %>%
  filter(date <= first_pos) %>%
  select(-first_pos) %>%
  distinct()


dat_for_age_model <- neg_line_list %>%
  filter(test_result == "positive") %>%
  filter(date >= "2020-03-14",
         date <= max(neg_line_list[["date"]]) - 14) %>%
  select(date, age) %>%
  mutate(date = as.numeric(date))

rcs_model <- ols(age ~ rcs(date,
                      quantile(date, c(0, .05, .25, .5, .75, .95, 1),
                               include.lowest = TRUE)),
            data = dat_for_age_model)


final_dat <- neg_line_list_filtered %>%
  count(date, test_result) %>%
  pivot_wider(names_from = test_result, values_from = n) %>%
  full_join(new_deaths_tbl) %>%
  right_join(., tibble(date = seq(min(.[["date"]]), max(.[["date"]]), by = 1))) %>%
  arrange(date) %>%
  replace(is.na(.), 0) %>%
  mutate(cases = positive,
         tests = negative + positive + other,
         lagged_pos_age = predict(rcs_model,
                                  newdata = tibble(date = as.numeric(date) - cases_age_lag))) %>%
  select(date, deaths, cases, tests, lagged_pos_age)

write_rds(final_dat, "data/oc_data.rds")
