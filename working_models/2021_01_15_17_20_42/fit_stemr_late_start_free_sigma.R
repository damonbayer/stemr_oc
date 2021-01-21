library(stemr)
library(extraDistr)
library(foreach)
library(doRNG)
library(doParallel)
library(tidyverse)
library(lubridate)
source("helper_functions.R")
# Load Data ---------------------------------------------------------------

if (Sys.info()[["sysname"]] == "Linux") {
  oc_data <- read_csv("/data/homezvol2/bayerd/uci_covid_modeling2/data/oc_data.csv")
} else if (Sys.info()[["sysname"]] == "Darwin") {
  oc_data <- read_csv("~/Documents/uci_covid_modeling2/data/oc_data.csv")
}

first_day <- ymd("2020-08-16")
last_day <- first_day + weeks(round(as.numeric((max(oc_data$date) - weeks(4) - first_day)) / 7)) - 1

dat <- oc_data %>%
  lump_oc_data(time_interval_in_days = 7,
               first_day,
               last_day)

obs_times <- dat$time

popsize <- 3175692L
log_popsize <- log(popsize)


# Build Model -------------------------------------------------------------
set.seed(12511)
strata <- NULL
compartments <- c("S", "E", "I", "R", "D")

rates <-
  list(
    rate(
      rate = "beta_t * I",
      from = "S",
      to = "E",
      incidence = F
    ),
    rate(
      rate = "gamma",
      from = "E",
      to = "I",
      incidence = T
    ),
    rate(
      rate = "mu_rec",
      from = "I",
      to = "R",
      incidence = F
    ),
    rate(
      rate = "mu_death",
      from = "I",
      to = "D",
      incidence = T
    )
  )

# n <- 1e5
# target_date <- ymd("2020-08-16")
# target_E_I <- oc_data %>% filter(date >= target_date - weeks(4), date <= target_date) %>% mutate(cases = cases * 7) %>% pull(cases) %>% sum()
# target_E <- round(target_E_I / 4)
# target_I <- target_E_I - target_E
# target_R <- round(0.115 * popsize)
# target_D <- oc_data %>% filter(date <= target_date) %>% pull(deaths) %>% sum()
# target_S <- popsize - (target_E + target_I + target_R + target_D)
# init_infected <- c(S = target_S, E = target_E, I = target_I, R = target_R, D = target_D)

init_states <- c(S = 2739459, E = 17488, I = 52463, R = 365205, D = 1077)

# my_function <- function(C) norm(as.matrix(quantile(x = extraDistr::rdirmnom(n = n, size = popsize, alpha = init_states / C)[,4], probs = c(0.025, 0.975)) - popsize * c(0.105, 0.124)), type = "F")
# res <- optimize(f = my_function, interval = c(100, 1000))
# res$minimum

C <- 732.389441820891

state_initializer <-
  list(
    stem_initializer(
      init_states = init_states,
      fixed = F,
      prior = init_states / C,
      dist = "dirmultinom"
    )
  )

adjacency <- NULL

tcovar <- dat %>% select(time, tests)

constants <- c(t0 = 0)
tmax <- max(obs_times)

foi_rw1 <- function(parameters, draws, log_pop = log_popsize) {
  log_R0_t <- numeric(length = length(draws))
  log_R0_t[1] <- log(parameters["R0_init"])

  for (t in 2:(length(log_R0_t))) {
    log_R0_t[t] <- log_R0_t[t - 1] + draws[t - 1] * parameters[["sigma"]]
  }

  exp(log_R0_t - log_pop + log(parameters[["mu_rec"]] + parameters[["mu_death"]]))
}


tparam <-
  list(tpar(
    tparam_name = "beta_t",
    times = c(0, head(obs_times,-1)),
    n_draws = length(obs_times),
    draws2par = foi_rw1
  ))

parameters_initial <-
  c(
    R0_init = 1,
    gamma = 1,
    mu_rec = 0.94,
    mu_death = 0.06,
    rho_death = 0.7,
    phi_death = 2.2,
    alpha0 = 4,
    alpha1 = 0.8,
    kappa = 2.2,
    sigma = 0.05
  )

dynamics <-
  stem_dynamics(
    rates = rates,
    tmax = tmax,
    parameters = parameters_initial,
    state_initializer = state_initializer,
    compartments = compartments,
    constants = constants,
    tcovar = tcovar,
    tparam = tparam,
    messages = T,
    compile_ode = T,
    compile_rates = F,
    compile_lna = F
  )

emissions <-
  list(
    emission(
      meas_var = "cases",
      distribution = "betabinomial",
      emission_params =
        c(
          "tests",
          "kappa * (exp(alpha0) * (E2I / popsize) ^ alpha1) / (exp(alpha0) * (E2I / popsize) ^ alpha1 + (1 - E2I / popsize) ^ alpha1)",
          "kappa * ((1 - E2I / popsize) ^ alpha1) / (exp(alpha0) * (E2I / popsize) ^ alpha1 + (1 - E2I / popsize) ^ alpha1)"
        ),
      incidence = TRUE,
      obstimes = obs_times
    ),
    emission(
      meas_var = "deaths",
      distribution = "negbinomial",
      emission_params = c("phi_death",
                          "rho_death * I2D"),
      incidence = TRUE,
      obstimes = obs_times
    ))

measurement_process <- stem_measure(emissions = emissions,
                                    dynamics = dynamics,
                                    data = select(dat, -ends_with("date"), -tests))

stem_object <-
  make_stem(dynamics = dynamics, measurement_process = measurement_process)

to_estimation_scale = function(params_nat) {
  c(R0_init_est = log(params_nat[["R0_init"]]), # log(R0_init)
    dur_latent_est = log(params_nat[["gamma"]]), # -log(dur_latent)
    dur_infec_est = log(params_nat[["mu_rec"]] + params_nat[["mu_death"]]), # -log(dur_infec)
    ifr_est = log(params_nat[["mu_death"]]) - log(params_nat[["mu_rec"]]), # logit(ifr)
    rho_death_est = logit(params_nat[["rho_death"]]), # logit(rho_death)
    phi_death_est = -0.5 * log(params_nat[["phi_death"]]), # -0.5 * log(phi_death)
    alpha0_est = log(params_nat[["alpha0"]]), # log(alpha0)
    alpha1_est = logit(params_nat[["alpha1"]]), # logit(alpha1)
    kappa_est = -0.5 * log(params_nat[["kappa"]]), # -0.5 * log(kappa)
    sigma_est = log(params_nat[["sigma"]])) # log(sigma)
}

# Need to account for sigma
from_estimation_scale = function(params_est) {
  c(R0_init = exp(params_est[["R0_init_est"]]),
    gamma = exp(params_est[["dur_latent_est"]]),
    mu_rec = exp(params_est[["dur_infec_est"]]) / (1 + exp(params_est[["ifr_est"]])),
    mu_death = exp(params_est[["dur_infec_est"]]) / (1 + exp(-params_est[["ifr_est"]])),
    rho_death = expit(params_est[["rho_death_est"]]),
    phi_death = exp(-2 * params_est[["phi_death_est"]]),
    alpha0 = exp(params_est[["alpha0_est"]]),
    alpha1 = expit(params_est[["alpha1_est"]]),
    kappa = exp(-2 * params_est[["kappa_est"]]),
    sigma = exp(params_est[["sigma_est"]]))
}

# Need to account for sigma
logprior =
  function(params_est) {
    sum(dnorm(params_est["R0_init_est"], -0.2554128198465173693599, 0.7, log = TRUE), # log(R0)
        dnorm(-params_est["dur_latent_est"],  log(3/7), 0.25, log = TRUE), # -log(dur_latent)
        dnorm(-params_est["dur_infec_est"], log(3), 0.175, log = TRUE), # -log(dur_infec)
        dbeta(expit(params_est["ifr_est"]), 2, 460, log = TRUE) + params_est["ifr_est"] - 2 * log(exp(params_est["ifr_est"]) + 1), # logit(ifr)
        dbeta(expit(params_est["rho_death_est"]), 8, 2, log = TRUE) + params_est["rho_death_est"] - 2 * log(exp(params_est["rho_death_est"]) + 1) , # logit(rho_death)
        dexp(exp(params_est["phi_death_est"]), 1, log = TRUE) +  params_est["phi_death_est"], # -0.5 * log(phi_death)
        dtnorm(exp(params_est["alpha0_est"]), mean = 4, sd = 2, a = 0, log = TRUE) + params_est["alpha0_est"], # log(alpha0)
        dbeta(expit(params_est["alpha1_est"]), 3, 1, log = TRUE) + params_est["alpha1_est"] - 2 * log(exp(params_est["alpha1_est"]) + 1), # logit(alpha1)
        dexp(exp(params_est["kappa_est"]), 1, log = T) +  params_est["kappa_est"], # -0.5 * log(kappa)
        dnorm(params_est["sigma_est"], -2.99573227714994, 0.421403565746124, log = T))
  }

priors <-
  list(logprior = logprior,
       to_estimation_scale = to_estimation_scale,
       from_estimation_scale = from_estimation_scale)

n_params <- length(parameters_initial)

par_initializer = function() {
  priors$from_estimation_scale(priors$to_estimation_scale(parameters_initial) +
                                 rnorm(n_params, 0, 0.1))
}

mcmc_kern <-
  mcmc_kernel(
    parameter_blocks =
      list(parblock(
        pars_nat = c("R0_init", "gamma", "mu_rec", "mu_death", "rho_death", "phi_death", "alpha0", "alpha1", "kappa", "sigma"),
        # par_est should have different names from pars_nat
        pars_est = c("R0_init_est", "dur_latent_est", "dur_infec_est", "ifr_est", "rho_death_est", "phi_death_est", "alpha0_est", "alpha1_est", "kappa_est", "sigma_est"),
        priors = priors,
        alg = "mvnss",
        sigma = diag(0.01, n_params),
        initializer = par_initializer,
        control =
          mvnss_control(stop_adaptation = 50000,
                        scale_cooling = 0.85,
                        scale_constant = 1,
                        step_size = 0.25))),
    lna_ess_control = lna_control(bracket_update_iter = 5e3,
                                  joint_initdist_update = FALSE),
    tparam_ess_control = tpar_control(bracket_update_iter = 1e4))

registerDoParallel(cores = future::availableCores())

n_chains <- 4
thinning_interval <- 100
iterations <- 350000

res <- foreach(chain = 1:n_chains,
               .packages = "stemr",
               .export = ls()) %dorng% {
                 fit_stem(stem_object = stem_object,
                          method = "ode", # or "lna"
                          mcmc_kern = mcmc_kern,
                          iterations = iterations,
                          thinning_interval = thinning_interval,
                          print_progress = 1e3)
               }

multi_chain_stem_fit <- list()
multi_chain_stem_fit$n_iterations <- iterations
multi_chain_stem_fit$n_chains <- n_chains
multi_chain_stem_fit$thinning_interval <- thinning_interval
multi_chain_stem_fit$stem_fit_list <- res
multi_chain_stem_fit$data <- dat
write_rds(multi_chain_stem_fit, paste0("res_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".rds"))
