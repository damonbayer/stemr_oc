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
# 5.4007449600867

deaths_od_sd <- sd(log(shape_draws))
dput(deaths_od_sd)
# 0.464859409320095

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

cases_model
pairs(cases_model)
phi_draws <- gather_draws(cases_model, phi)$.value


cases_od_mean <- mean(log(phi_draws))
# dput(cases_od_mean)
6.99157747978605

cases_od_sd <- sd(log(phi_draws))
# dput(cases_od_sd)
0.237796375571439

exp(qnorm(p = desired_quantiles, mean = cases_od_mean, cases_od_sd))
quantile(phi_draws, desired_quantiles)
