library(tidyverse)
library(tidybayes)
library(scales)
library(cowplot)
library(fs)
library(stemr)
library(coda)
library(extraDistr)
source('stemr_functions.R')
theme_set(cowplot::theme_minimal_grid())

model_names <-
  dir_ls("/Users/damon/Documents/stemr_oc/final_models/results/", recurse = F) %>%
  path_ext_remove() %>%
  path_file() %>%
  enframe("model_name", name = NULL) %>%
  filter(str_starts(model_name, "20", negate = T)) %>%
  deframe()

for (model_name in model_names) {

print(model_name)
multi_chain_stem_fit <- read_rds(path("final_models", "results", model_name, ext = "rds"))

popsize <- multi_chain_stem_fit$stem_fit_list[[1]]$dynamics$popsize
log_popsize <- log(popsize)

time_date_conversion <- select(multi_chain_stem_fit$data, time, date = end_date) %>%
  add_row(., time = 0, date = min(.$date) - 7, .before = T)

# Gather Parameters -------------------------------------------------------
epi_curves_p <- extract_epi_curves(multi_chain_stem_fit = multi_chain_stem_fit, curve_type = "p", tidy = F)

stemr_tparams <-
  imap_dfr(multi_chain_stem_fit$stem_fit_list, ~{
    cbind(
      gather_array(.x$results$posterior$tparam_samples,
                   value = value, time, var, .iteration),
      .chain = .y)}) %>%
  as_tibble() %>%
  left_join(enframe(c("beta_t", "ifr_t", "alpha0_t"), name = "var", value = "name"), by = "var") %>%
  mutate(time = time - 1) %>%
  select(-var) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  mutate(.draw = tidybayes:::draw_from_chain_and_iteration_(.chain, .iteration)) %>%
  select(.chain, .iteration, .draw, time, everything())

stemr_transitions <-
  imap_dfr(multi_chain_stem_fit$stem_fit_list,  ~{
    .x$results$posterior$latent_paths %>%
      gather_array(value, row, name, .iteration) %>%
      as_tibble() %>%
      mutate(name = dimnames(.x$results$posterior$latent_paths)[[2]][name],
             time = .x$results$posterior$latent_paths[,1,1][row]) %>%
      filter(name != "time") %>%
      select(.iteration, time, name, value, -row) %>%
      pivot_wider(names_from = name, values_from = value) %>%
      mutate(.chain = .y)
  }) %>%
  mutate(.draw = tidybayes:::draw_from_chain_and_iteration_(.chain, .iteration)) %>%
  select(.chain, .iteration, .draw, time, everything())

stemr_all_params <-
  left_join(stemr_tparams, stemr_transitions,
            by = c(".chain", ".iteration", ".draw", "time")) %>%
  left_join(pivot_wider(extract_epi_curves(multi_chain_stem_fit = multi_chain_stem_fit, curve_type = "p", tidy = T), values_from = value, names_from = name),
            by = c(".chain", ".iteration", ".draw", "time")) %>%
  left_join(tidy_draws(
    extract_stem_parameter_posterior(multi_chain_stem_fit = multi_chain_stem_fit,
                                     transform = "natural")),
    by = c(".chain", ".iteration", ".draw")) %>%
  left_join(multi_chain_stem_fit$data,
            by = "time") %>%
  left_join(time_date_conversion, by = "time")

rm(stemr_tparams, stemr_transitions)

# PP Data -----------------------------------------------------------------
set.seed(1)
pp_data <-
  stemr_all_params %>%
  filter(time != 0) %>%
  mutate(tests = tests,
         alpha = kappa * (exp(alpha0_t) * (E2I / popsize)) / (exp(alpha0_t) * (E2I / popsize) + (1 - E2I / popsize)),
         beta = kappa * ((1 - E2I / popsize)) / ((E2I / popsize) + (1 - E2I / popsize)),
         size = phi_death,
         mu = rho_death * I2D,
         seroprev_tests = seroprev_tests,
         seroprev_prob = R / (S + E + I + R)) %>%
  mutate(
    cases = pmap_dbl(
      list(tests = as.list(tests),
           alpha = as.list(alpha),
           beta = as.list(beta)),
      function(tests, alpha, beta) rbbinom(n = 1, size = tests, alpha = alpha, beta = beta)),
    deaths = pmap_dbl(
      list(mu = as.list(mu),
           size = as.list(size)),
      function(mu, size) rnbinom(n = 1, size = size, mu = mu)),
    # function(mu, size) rnbinom(n = 1, size = size, mu = mu))
    seroprev_cases = pmap_dbl(
      list(seroprev_tests = as.list(seroprev_tests),
           seroprev_prob = as.list(seroprev_prob)),
      function(seroprev_tests, seroprev_prob) ifelse(is.na(seroprev_tests), NA, rbinom(n = 1, size = seroprev_tests, prob = seroprev_prob)))
  ) %>%
  select(starts_with("."), date = end_date, tests, cases, deaths, seroprev_tests, seroprev_cases) %>%
  mutate(pos = cases / tests,
         seroprev_pos = seroprev_cases / seroprev_tests) %>%
  select(starts_with("."), date, deaths, pos, seroprev_pos) %>%
  pivot_longer(c(deaths, pos, seroprev_pos)) %>%
  drop_na()

# Save Data ---------------------------------------------------------------

write_rds(stemr_all_params, path("final_models", "results", "stemr_all_params", str_c(model_name, "_stemr_all_params"), ext = "rds"))
write_rds(pp_data, path("final_models", "results", "pp_data", str_c(model_name, "_pp_data"), ext = "rds"))

}
