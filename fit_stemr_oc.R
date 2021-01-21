library(stemr)
library(extraDistr)
library(foreach)
library(doRNG)
library(doParallel)
library(tidyverse)

# Load Data ---------------------------------------------------------------
popsize <- 3068927L
log_popsize <- log(popsize)

oc_data <- read_rds("data/oc_data_no_snf.rds")

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
      dplyr::summarize(start_date = min(date),
                       end_date = max(date),
                       cases = sum(cases),
                       tests = sum(tests),
                       deaths = sum(deaths)) %>%
      dplyr::select(-lump) %>%
      arrange(start_date)
  }

time_interval_in_days <- 7

lumped_oc_data <- lump_oc_data(oc_data, time_interval_in_days = time_interval_in_days, first_day = "2020-03-30", last_day = "2020-10-11") %>%
  left_join(select(oc_data, end_date = date, lagged_pos_ifr))

dat <- lumped_oc_data %>%
  mutate(time = as.numeric((lumped_oc_data$end_date - min(lumped_oc_data$start_date) + 1L) / time_interval_in_days)) %>%
  select(time, cases, tests, deaths, lagged_pos_ifr) %>%
  full_join(tibble(time = as.numeric((ymd("2020-8-16") - min(lumped_oc_data$start_date) + 1L) / time_interval_in_days),
                   prevalence = as.integer(0.115 * popsize))) %>%
  arrange(time)

obs_times <- dat$time

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
      rate = "(1 - ifr_t) / dur_infec",
      from = "I",
      to = "R",
      incidence = F
    ),
    rate(
      rate = "ifr_t / dur_infec",
      from = "I",
      to = "D",
      incidence = T
    )
  )

init_infected <-
  sum(oc_data$cases[oc_data$date >= "2020-03-15" & oc_data$date <= "2020-03-30"]) * 7

init_states <-
  c(S = popsize - init_infected,
    E = round(init_infected / 3),
    I = init_infected - round(init_infected / 3),
    R = 0,
    D = 0)

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

tcovar <- dat %>% mutate(logit_lagged_pos_ifr = qlogis(lagged_pos_ifr)) %>% select(time, tests, logit_lagged_pos_ifr)

constants <- c(t0 = 0, sigma = 0.05)
tmax <- max(obs_times)

foi_rw1 <- function(parameters, draws, log_pop = log_popsize) {
  log_R0_t <- numeric(length = length(draws))
  log_R0_t[1] <- log(parameters["R0_init"])

  for (t in 2:(length(log_R0_t))) {
    log_R0_t[t] <- log_R0_t[t - 1] + draws[t - 1] * constants[["sigma"]]
  }

  exp(log_R0_t - log_pop - log(1 / parameters[["dur_infec"]]))
}

ifr_fn <- function(parameters, draws, ifr_est = unname(tcovar[["logit_lagged_pos_ifr"]])) {
  plogis(ifr_est + parameters[["ifr_shift"]])
}


tparam <-
  list(tpar(
    tparam_name = "beta_t",
    times = c(0, head(obs_times,-1)),
    n_draws = length(obs_times),
    draws2par = foi_rw1
  ),
  tpar(tparam_name = "ifr_t",
       times = c(0, head(obs_times,-1)),
       n_draws = 1, # would like to make this 0 at some point
       draws2par = ifr_fn))

parameters_initial <-
  c(
    R0_init = 2,
    gamma = 1,
    dur_infec = 0.94,
    ifr_shift = 0.06,
    rho_death = 0.7,
    phi_death = 2.2,
    alpha0 = 4,
    alpha1 = 0.8,
    kappa = 2.2
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

prevalence_emission_sd <- optimize(f = function(x) norm(as.matrix(qnorm(p = c(0.025, 0.5, 0.975), mean = popsize * 0.115, sd = x) - popsize * c(0.105, 0.115, 0.124)), type = "f"), lower = 10000, upper = 30000)$minimum

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
    ),
    emission(meas_var = "prevalence",
             distribution = "gaussian",
             emission_params = c("R", as.character(prevalence_emission_sd)),
             incidence = F,
             obstimes = obs_times
  ))

measurement_process <- stem_measure(emissions = emissions,
                                    dynamics = dynamics,
                                    data = select(dat, -tests, -lagged_pos_ifr))

stem_object <-
  make_stem(dynamics = dynamics, measurement_process = measurement_process)

# ifr_shift and dur_infec
to_estimation_scale = function(params_nat) {
  c(R0_init_est = log(params_nat[["R0_init"]]), # log(R0_init)
    dur_latent_est = log(params_nat[["gamma"]]), # -log(dur_latent)
    dur_infec_est = log(params_nat[["dur_infec"]]), # log(dur_infec)
    ifr_shift_est = params_nat[["ifr_shift"]], # ifr_shift
    rho_death_est = logit(params_nat[["rho_death"]]), # logit(rho_death)
    phi_death_est = -0.5 * log(params_nat[["phi_death"]]), # -0.5 * log(phi_death)
    alpha0_est = log(params_nat[["alpha0"]]), # log(alpha0)
    alpha1_est = logit(params_nat[["alpha1"]]), # logit(alpha1)
    kappa_est = -0.5 * log(params_nat[["kappa"]])) # -0.5 * log(kappa)
}

# Need to account for sigma
from_estimation_scale = function(params_est) {
  c(R0_init = exp(params_est[["R0_init_est"]]),
    gamma = exp(params_est[["dur_latent_est"]]),
    dur_infec = exp(params_est[["dur_infec_est"]]),
    ifr_shift = params_est[["ifr_shift_est"]],
    rho_death = expit(params_est[["rho_death_est"]]),
    phi_death = exp(-2 * params_est[["phi_death_est"]]),
    alpha0 = exp(params_est[["alpha0_est"]]),
    alpha1 = expit(params_est[["alpha1_est"]]),
    kappa = exp(-2 * params_est[["kappa_est"]]))
}

# Need to account for sigma
logprior =
  function(params_est) {
    sum(dnorm(params_est["R0_init_est"], -0.2554128198465173693599, 0.7, log = TRUE), # log(R0)
        dnorm(-params_est["dur_latent_est"],  log(3/7), 0.25, log = TRUE), # -log(dur_latent)
        dnorm(-params_est["dur_infec_est"], log(3), 0.175, log = TRUE), # -log(dur_infec)
        dnorm(params_est["ifr_shift_est"], 0, 0.4, log = TRUE), # logit(ifr)
        dbeta(expit(params_est["rho_death_est"]), 8, 2, log = TRUE) + params_est["rho_death_est"] - 2 * log(exp(params_est["rho_death_est"]) + 1) , # logit(rho_death)
        dexp(exp(params_est["phi_death_est"]), 1, log = TRUE) +  params_est["phi_death_est"], # -0.5 * log(phi_death)
        dtnorm(exp(params_est["alpha0_est"]), mean = 4, sd = 2, a = 0, log = TRUE) + params_est["alpha0_est"], # log(alpha0)
        dbeta(expit(params_est["alpha1_est"]), 3, 1, log = TRUE) + params_est["alpha1_est"] - 2 * log(exp(params_est["alpha1_est"]) + 1), # logit(alpha1)
        dexp(exp(params_est["kappa_est"]), 1, log = T) +  params_est["kappa_est"]) # -0.5 * log(kappa)
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
        pars_nat = c("R0_init", "gamma", "dur_infec", "ifr_shift", "rho_death", "phi_death", "alpha0", "alpha1", "kappa"),
        pars_est = c("R0_init_est", "dur_latent_est", "dur_infec_est", "ifr_shift_est", "rho_death_est", "phi_death_est", "alpha0_est", "alpha1_est", "kappa_est"),
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
