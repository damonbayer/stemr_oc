library(stemr)
library(extraDistr)
library(foreach)
library(doRNG)
library(doParallel)
library(tidyverse)

# Load Data ---------------------------------------------------------------
oc_data <-
  read_csv(
    "/Users/damon/Documents/uci_covid_modeling/data/oc_data.csv",
    # "oc_data.csv",
    col_types = cols(
      X1 = col_skip(),
      posted_date = col_date(format = "%Y-%m-%d"),
      new_cases = col_integer(),
      new_tests = col_integer(),
      new_deaths = col_integer(),
      new_deaths_calcat = col_skip()
    )
  ) %>%
  rename(date = posted_date)

lump_oc_data <-
  function(oc_data,
           time_interval_in_days,
           first_day = "0000-01-01",
           last_day = "9999-12-31") {
    oc_data %>%
    filter(date >= lubridate::ymd(first_day),
           date <= lubridate::ymd(last_day)) %>%
    group_by(lump = as.integer(floor((max(date) - date) / time_interval_in_days))) %>%
    filter(n() == time_interval_in_days) %>%
    summarize(start_date = min(date),
              end_date = max(date),
              new_cases = sum(new_cases),
              new_tests = sum(new_tests),
              new_deaths = sum(new_deaths)) %>%
    dplyr::select(-lump) %>%
    arrange(start_date)
}

lumped_oc_data <- lump_oc_data(oc_data, time_interval_in_days = 7, first_day = "2020-03-15", last_day = "2020-10-15")

dat <- lumped_oc_data %>%
  mutate(time = as.numeric((end_date - min(start_date) + 1L) / 7)) %>%
  select(time, cases = new_cases, deaths = new_deaths, tests = new_tests)

obs_times <- dat$time

popsize <- 3068927L
log_popsize <- log(popsize)


# Build Model -------------------------------------------------------------
set.seed(12511)
strata <- NULL
compartments <- c("S", "E", "Ie", "Ip", "R", "D")

rates <-
  list(
    rate(
      rate = "beta_t * (Ie + 0.8 * Ip)",
      from = "S",
      to = "E",
      incidence = F
    ),
    rate(
      rate = "gamma",
      from = "E",
      to = "Ie",
      incidence = F
    ),
    rate(
      rate = "nu_early",
      from = "Ie",
      to = "Ip",
      incidence = T
    ),
    rate(
      rate = "mu_rec",
      from = "Ip",
      to = "R",
      incidence = F
    ),
    rate(
      rate = "mu_death",
      from = "Ip",
      to = "D",
      incidence = T
    )
  )

init_infected <-
  sum(oc_data$new_cases[oc_data$date <= "2020-03-14"]) * 10

init_states <-
  c(popsize - init_infected, lengths(parallel::splitIndices(init_infected, 3)), 0, 0) %>%
  `names<-`(compartments)

state_initializer <-
  list(
    stem_initializer(
      init_states = init_states,
      fixed = F,
      prior = init_states / 10,
      dist = "dirmultinom"
    )
  )

adjacency <- NULL

tcovar <- dat %>% select(time, tests)

parameters <-
  c(
    log_R0_init = log(2),
    gamma = 1,
    nu_early = 1,
    mu_rec = 0.94,
    mu_death = 0.06,
    rho_death = 0.7,
    phi_death = 2.2,
    alpha0 = 4,
    alpha1 = 0.8,
    kappa = 2.2,
    sigma = 0.05
  )

constants <- c(t0 = 0)
tmax <- max(obs_times)

foi_rw1 <- function(parameters, draws, log_pop = log_popsize) {
  log_R0_t <- numeric(length = length(draws))
  log_R0_t[1] <- parameters["log_R0_init"]

  for (t in 2:(length(log_R0_t))) {
    log_R0_t[t] <- log_R0_t[t - 1] + draws[t - 1] * parameters["sigma"]
  }

  exp(log_R0_t - log_pop - log(1 / parameters[["nu_early"]] + 0.8 / (parameters[["mu_rec"]] + parameters[["mu_death"]])))
}

tparam <-
  list(tpar(
    tparam_name = "beta_t",
    times = c(0, head(obs_times,-1)),
    n_draws = length(obs_times),
    draws2par = foi_rw1
  ))

dynamics <-
  stem_dynamics(
    rates = rates,
    tmax = tmax,
    parameters = parameters,
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
          "kappa * (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1) / (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1 + (1 - Ie2Ip/popsize) ^ alpha1)",
          "kappa * ((1 - Ie2Ip/popsize) ^ alpha1) / (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1 + (1 - Ie2Ip/popsize) ^ alpha1)"
        ),
      incidence = TRUE,
      obstimes = obs_times
    ),
    emission(
      meas_var = "deaths",
      distribution = "negbinomial",
      emission_params = c("phi_death",
                          "rho_death * Ip2D"),
      incidence = TRUE,
      obstimes = obs_times
    )
  )

measurement_process <- stem_measure(emissions = emissions,
                                    dynamics = dynamics,
                                    data = select(dat, -tests))

stem_object <-
  make_stem(dynamics = dynamics, measurement_process = measurement_process)


parameters <-
  c(
    log_R0_init = log(2),
    gamma = 1,
    nu_early = 1,
    mu_rec = 0.94,
    mu_death = 0.06,
    rho_death = 0.7,
    phi_death = 2.2,
    alpha0 = 4,
    alpha1 = 0.8,
    kappa = 2.2,
    sigma = 0.05
  )

# Need to account for sigma
to_estimation_scale = function(params_nat) {
  c(log_R0_init_est = params_nat[["log_R0_init"]],
    dur_latent_est = log(params_nat[["gamma"]]), # -log(dur_latent)
    dur_early_est = log(params_nat[["nu_early"]]), # -log(dur_early)
    dur_progress_est = log(params_nat[["mu_rec"]] + params_nat[["mu_death"]]), # -log(dur_progress)
    ifr_est = log(params_nat[["mu_death"]]) - log(params_nat[["mu_rec"]]), # logit(ifr)
    rho_death_est = logit(params_nat[["rho_death"]]), # logit(rho_death)
    phi_death_est = -0.5 * log(params_nat[["phi_death"]]), # -0.5 * log(phi_death)
    alpha0_est = log(params_nat[["alpha0"]]), # log(alpha0)
    alpha1_est = logit(params_nat[["alpha1"]]), # logit(alpha1)
    kappa_est = -0.5 * log(params_nat[["kappa"]]), # -0.5 * log(kappa)
    sigma_est = log(params_nat[["sigma"]]))
}

# Need to account for sigma
from_estimation_scale = function(params_est) {
  c(log_R0_init = params_est[["log_R0_init_est"]],
    gamma = exp(params_est[["dur_latent_est"]]),
    nu_early = exp(params_est[["dur_early_est"]]),
    mu_rec = exp(params_est[["dur_progress_est"]]) / (1 + exp(params_est[["ifr_est"]])),
    mu_death = exp(params_est[["dur_progress_est"]]) / (1 + exp(-params_est[["ifr_est"]])),
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
    sum(dnorm(params_est["log_R0_init_est"], -0.2554128198465173693599, 0.7, log = TRUE), # log(R0)
        dnorm(-params_est["dur_latent_est"], 0, 0.22, log = TRUE), # -log(dur_latent)
        dnorm(-params_est["dur_early_est"], 0, 0.22, log = TRUE), # -log(dur_early)
        dnorm(-params_est["dur_progress_est"], 0, 0.22, log = TRUE), # -log(dur_progress)
        dbeta(expit(params_est["ifr_est"]), 1.5, 200, log = TRUE) + params_est["ifr_est"] - 2 * log(exp(params_est["ifr_est"]) + 1), # logit(ifr)
        dbeta(expit(params_est["rho_death_est"]), 8, 2, log = TRUE) + params_est["rho_death_est"] - 2 * log(exp(params_est["rho_death_est"]) + 1) , # logit(rho_death)
        dexp(exp(params_est["phi_death_est"]), 1, log = TRUE) +  params_est["phi_death_est"], # -0.5 * log(phi_death)
        dtnorm(exp(params_est["alpha0_est"]), mean = 4, sd = 2, a = 0, log = TRUE) + params_est["alpha0_est"], # log(alpha0)
        dbeta(expit(params_est["alpha1_est"]), 3, 1, log = TRUE) + params_est["alpha1_est"] - 2 * log(exp(params_est["alpha1_est"]) + 1), # logit(alpha1)
        dexp(exp(params_est["kappa_est"]), 1, log = T) +  params_est["kappa_est"], # -0.5 * log(kappa)
        dexp(exp(params_est[["sigma_est"]]), 20, log = TRUE) + params_est[["sigma_est"]])
  }

priors <-
  list(logprior = logprior,
       to_estimation_scale = to_estimation_scale,
       from_estimation_scale = from_estimation_scale)

n_params <- length(parameters)

par_initializer = function() {
  priors$from_estimation_scale(priors$to_estimation_scale(parameters) +
                                 rnorm(n_params, 0, 0.1))
}

mcmc_kern <-
  mcmc_kernel(
    parameter_blocks =
      list(parblock(
        pars_nat = c("log_R0_init", "gamma", "nu_early", "mu_rec", "mu_death", "rho_death", "phi_death", "alpha0", "alpha1", "kappa", "sigma"),
        # par_est should have different names from pars_nat
        pars_est = c("log_R0_init_est", "dur_latent_est", "dur_early_est", "dur_progress_est", "ifr_est", "rho_death_est", "phi_death_est", "alpha0_est", "alpha1_est", "kappa_est", "sigma_est"),
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

res <- foreach(chain = 1:4,
               .packages = "stemr",
               .export = ls()) %dorng% {
                 fit_stem(stem_object = stem_object,
                          method = "ode", # or "lna"
                          mcmc_kern = mcmc_kern,
                          iterations = iterations,
                          thinning_interval = thinning_interval,
                          print_progress = 1e3)
               }


write_rds(res, paste0("res_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".rds"))
