library(tidyverse)
library(brms)
library(tidybayes)
library(bayesplot)

options(mc.cores = parallel::detectCores())
dat <- read_csv("final_models/data/simulated_data.csv")

desired_quantiles <- c(0.05, 0.5, 0.95)

deaths_model <- brm(bf(deaths ~ s(time)),
             data = dat,
             family = negbinomial(),
             cores = 4,
             seed = 17,
             iter = 4000,
             warmup = 1000,
             # thin = 10,
             # refresh = 0
             # control = list(adapt_delta = 0.99),
             control = list(adapt_delta = 0.99, max_treedepth = 12),
             backend = "cmdstanr")



shape_draws <- gather_draws(deaths_model, shape)$.value
deaths_od_mean <- mean(log(shape_draws))
dput(deaths_od_mean)
deaths_od_mean <- 5.4007449600867

deaths_od_sd <- sd(log(shape_draws))
dput(deaths_od_sd)
deaths_od_sd <- 0.464859409320095

exp(qnorm(p = desired_quantiles, mean = deaths_od_mean, deaths_od_sd))
quantile(shape_draws, desired_quantiles)

# https://cran.r-project.org/web/packages/brms/vignettes/brms_customfamilies.html

beta_binomial2 <- custom_family(
  "beta_binomial2", dpars = c("mu", "phi"),
  links = c("logit", "log"), lb = c(NA, 0),
  type = "int", vars = "vint1[n]"
)

stan_funs <- "
  real beta_binomial2_lpmf(int y, real mu, real phi, int T) {
    return beta_binomial_lpmf(y | T, mu * phi, (1 - mu) * phi);
  }
  int beta_binomial2_rng(real mu, real phi, int T) {
    return beta_binomial_rng(T, mu * phi, (1 - mu) * phi);
  }
"

stanvars <- stanvar(scode = stan_funs, block = "functions")


cases_model <- brm(
  cases | vint(tests) ~ s(time), data = dat, 
  family = beta_binomial2, stanvars = stanvars,
  cores = 4,
  seed = 18,
  iter = 8000,
  warmup = 4000,
  # control = list(adapt_delta = 0.8)
)

pairs(cases_model)
phi_draws <- gather_draws(cases_model, phi)$.value


cases_od_mean <- mean(log(phi_draws))
# dput(cases_od_mean)
cases_od_mean <- 6.99157747978605

cases_od_sd <- sd(log(phi_draws))
# dput(cases_od_sd)
cases_od_sd <- 0.237796375571439

exp(qnorm(p = desired_quantiles, mean = cases_od_mean, cases_od_sd))
quantile(phi_draws, desired_quantiles)


library(extraDistr)
stemr_results <- read_rds("final_models/results/stemr_all_params/2021-02-28_stemr_all_params.rds")

values_for_comparison <- 
  bind_rows(
    stemr_results %>% 
      filter(time == 0) %>% 
      select(phi_death, kappa) %>% 
      rename(phi_cases = kappa) %>% 
      pivot_longer(everything()) %>% 
      mutate(model = "stemr"),
    bind_rows(
      enframe(shape_draws, name = NULL) %>% 
        mutate(name = "phi_death"),
      enframe(phi_draws, name = NULL) %>% 
        mutate(name = "phi_cases")) %>% 
      mutate(model = "brms"))


values_for_comparison %>% 
  filter(name == "phi_death") %>% 
  ggplot(aes(value, group = model, fill = model)) +
  # facet_wrap(. ~ name, scales = "free") +
  stat_halfeye(normalize = "panels", alpha = 0.5) +
  cowplot::theme_cowplot()

values_for_comparison %>% 
  filter(name == "phi_cases") %>% 
  ggplot(aes(value, group = model, fill = model)) +
  # facet_wrap(. ~ name, scales = "free") +
  stat_halfeye(normalize = "panels", alpha = 0.5) +
  cowplot::theme_cowplot()

# values_for_comparison %>% 
# filter(name == "phi_cases") %>% 
#   group_by(model) %>% 
#   summarize(mean(value), min(value), max(value))


rbeta_binomial2 <- function(n, size, mu, phi) {
  alpha <- mu * phi
  beta <- (1 - mu) * phi
  rbbinom(n = n, size = size, alpha = alpha, beta = beta)
}

post_pred_cases <- 
  values_for_comparison %>%
  filter(name == "phi_cases") %>% 
  crossing(tests = c(5000, 10000, 50000, 100000, 150000),
           positivity = c(0.01, 0.05, 0.1, 0.2)) %>%
  mutate(positive_tests = pmap_dbl(.l = list(size = tests,
                                             mu = positivity,
                                             phi = value),
                                   function(size, mu, phi)
                                     rbeta_binomial2(n = 1, size = size, mu = mu, phi = phi)))

post_pred_cases %>% 
  mutate(observed_positivity = positive_tests / tests) %>% 
  ggplot(aes(observed_positivity, group = model, fill = model)) +
  facet_grid(tests ~ positivity, scales = "free", labeller = label_both) +
  stat_halfeye(normalize = "panels", alpha = 0.5) +
  cowplot::theme_cowplot() +
  scale_x_continuous(labels = ~scales::percent(., accuracy = 1))


post_pred_deaths <-
  values_for_comparison %>%
  filter(name == "phi_death") %>% 
  crossing(mean_deaths = c(10, 50, 100, 200, 400)) %>% 
  mutate(observed_deaths = map2_dbl(.x = mean_deaths, .y = value, ~rnbinom(n = 1, size = .y, mu = .x)))



post_pred_deaths %>% 
  ggplot(aes(observed_deaths, group = model, fill = model)) +
  facet_wrap(. ~ mean_deaths, scales = "free", labeller = label_both) +
  stat_halfeye(normalize = "panels", alpha = 0.5) +
  cowplot::theme_cowplot()
