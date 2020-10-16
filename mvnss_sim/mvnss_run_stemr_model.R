library(tidyverse)
library(stemr)
library(extraDistr)
library(foreach)
library(doRNG)
library(doParallel)
library(furrr)
registerDoParallel(cores = future::availableCores())
theme_set(cowplot::theme_cowplot())


# Read Data ---------------------------------------------------------------
dat <- read_rds("~/Documents/uci_covid_modeling/code/SEIeIpRD/oc/2020-08-12_2020-09-16/model_objects.rds")$lumped_ochca_covid

init_states <- c(S = 3012564L, E = 14866L, Ie = 20907L, Ip = 20590L, R = 0L, D = 0L)
popsize <- sum(init_states)

# Build stem_object -------------------------------------------------------
obs_times <- as.numeric((dat$end_date - min(dat$start_date) + 1L) / 7)

strata <- NULL # no strata
compartments <- c("S", "E", "Ie", "Ip", "R", "D")

rates <-
  list(rate(rate = "beta * (Ie + 0.8 * Ip)", # individual level rate (unlumped)
            from = "S",        # source compartment
            to   = "E",        # destination compartment
            incidence = F),    # compute incidence of S2I transitions, required for simulating incidence data
       rate(rate = "gamma",
            from = "E",
            to = "Ie",
            incidence = F),
       rate(rate = "nu_early",
            from = "Ie",
            to = "Ip",
            incidence = T),
       rate(rate = "mu_rec",
            from = "Ip",
            to = "R",
            incidence = F),
       rate(rate = "mu_death",       # individual level rate
            from = "Ip",        # source compartment
            to   = "D",        # destination compartment
            incidence = TRUE)) # compute incidence of I2R transitions (not required for simulating data)

state_initializer <-
  list(stem_initializer(
    init_states = init_states, # must match compartment names
    fixed = T,
    # prior = init_states / C,
    dist = "dirmultinom"
  )) # initial state fixed for simulation, we'll change this later

parameters <- numeric(10); names(parameters) <- c("beta", "gamma", "nu_early", "mu_rec", "mu_death", "rho_death", "phi_death", "alpha0", "alpha1", "kappa")
constants <- c(t0 = 0)
tcovar <- data.frame(time = obs_times,
                     tests = dat$new_tests)
tmax <- max(tcovar$time)

# list of emission distribution lists (analogous to rate specification)
emissions <-
  list(emission(meas_var = "cases", # transition or compartment being measured (S->I transitions)
                distribution    = "betabinomial",        # emission distribution
                emission_params =
                  c("tests",
                    "kappa * (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1) / (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1 + (1 - Ie2Ip/popsize) ^ alpha1)",
                    "kappa * ((1 - Ie2Ip/popsize) ^ alpha1) / (exp(alpha0) * (Ie2Ip / popsize) ^ alpha1 + (1 - Ie2Ip/popsize) ^ alpha1)"), # distribution pars, here overdispersion and mean
                incidence       = TRUE,                  # is the data incidence
                obstimes        = obs_times), # vector of observation times
       emission(meas_var = "deaths",
                distribution = "negbinomial",
                emission_params = c("phi_death",
                                    "rho_death * Ip2D"),
                incidence = T,
                obstimes        = obs_times)) # vector of observation times)
# list of emission distribution lists (analogous to rate specification)

dynamics <-
  stem_dynamics(
    rates = rates,
    parameters = parameters,
    state_initializer = state_initializer,
    compartments = compartments,
    constants = constants,
    tcovar = tcovar,
    tmax = tmax,
    compile_ode = T,   # compile ODE functions
    compile_rates = F, # compile MJP functions for Gillespie simulation
    compile_lna = T,   # compile LNA functions
    messages = F       # don't print messages
  )



measurement_process <-
  stem_measure(emissions = emissions,
               dynamics = dynamics,
               data = dat %>%
                 mutate(t= obs_times) %>%
                 select(t, cases = new_cases, deaths = new_deaths))

stem_object <- make_stem(dynamics = dynamics, measurement_process = measurement_process)



# Build Priors ------------------------------------------------------------
to_estimation_scale = function(params_nat) {
  c(R0_est = log(params_nat[["beta"]]) + log(popsize) + log(1 / params_nat[["nu_early"]] + 0.8 / (params_nat[["mu_rec"]] + params_nat[["mu_death"]])), # log(R0)
    dur_latent_est = log(params_nat[["gamma"]]), # -log(dur_latent)
    dur_early_est = log(params_nat[["nu_early"]]), # -log(dur_early)
    dur_progress_est = log(params_nat[["mu_rec"]] + params_nat[["mu_death"]]), # -log(dur_progress)
    ifr_est = log(params_nat[["mu_death"]]) - log(params_nat[["mu_rec"]]), # logit(ifr)
    rho_death_est = logit(params_nat[["rho_death"]]), # logit(rho_death)
    phi_death_est = -0.5 * log(params_nat[["phi_death"]]), # -0.5 * log(phi_death)
    alpha0_est = log(params_nat[["alpha0"]]), # log(alpha0)
    alpha1_est = logit(params_nat[["alpha1"]]), # logit(alpha1)
    kappa_est = -0.5 * log(params_nat[["kappa"]])) # -0.5 * log(kappa)
}


from_estimation_scale = function(params_est) {
  c(beta = exp(params_est[["R0_est"]] - log(popsize) - log(exp(-params_est[["dur_early_est"]]) + 0.8 * exp(-params_est[["dur_progress_est"]]))),
    gamma = exp(params_est[["dur_latent_est"]]),
    nu_early = exp(params_est[["dur_early_est"]]),
    mu_rec = exp(params_est[["dur_progress_est"]]) / (1 + exp(params_est[["ifr_est"]])),
    mu_death = exp(params_est[["dur_progress_est"]]) / (1 + exp(-params_est[["ifr_est"]])),
    rho_death = expit(params_est[["rho_death_est"]]),
    phi_death = exp(-2 * params_est[["phi_death_est"]]),
    alpha0 = exp(params_est[["alpha0_est"]]),
    alpha1 = expit(params_est[["alpha1_est"]]),
    kappa = exp(-2 * params_est[["kappa_est"]]))
}


logprior =
  function(params_est) {
    sum(dnorm(params_est["R0_est"], -0.25, 0.7, log = TRUE), # log(R0)
        dnorm(-params_est["dur_latent_est"], 0, 0.22, log = TRUE), # -log(dur_latent)
        dnorm(-params_est["dur_early_est"], 0, 0.22, log = TRUE), # -log(dur_early)
        dnorm(-params_est["dur_progress_est"], 0, 0.22, log = TRUE), # -log(dur_progress)
        dbeta(expit(params_est["ifr_est"]), 1.5, 200, log = TRUE) + params_est["ifr_est"] - 2 * log(exp(params_est["ifr_est"]) + 1), # logit(ifr)
        dbeta(expit(params_est["rho_death_est"]), 8, 2, log = TRUE) + params_est["rho_death_est"] - 2 * log(exp(params_est["rho_death_est"]) + 1) , # logit(rho_death)
        dexp(exp(params_est["phi_death_est"]), 1, log = TRUE) +  params_est["phi_death_est"], # -0.5 * log(phi_death)
        dtnorm(exp(params_est["alpha0_est"]), mean = 4, sd = 2, a = 0, log = TRUE) + params_est["alpha0_est"], # log(alpha0)
        dbeta(expit(params_est["alpha1_est"]), 3, 1, log = TRUE) + params_est["alpha1_est"] - 2 * log(exp(params_est["alpha1_est"]) + 1), # logit(alpha1)
        dexp(exp(params_est["kappa_est"]), 1, log = T) +  params_est["kappa_est"]) # -0.5 * log(kappa)
  }

priors <- list(logprior = logprior,
               to_estimation_scale = to_estimation_scale,
               from_estimation_scale = from_estimation_scale)


# Build par_initializer ---------------------------------------------------
true_pars =
  c(R0       = 1.5,    # basic reproduction number
    dur_latent = 1,
    dur_early   = 1,      # infectious period duration = 2 days
    dur_progress = 1,
    ifr = 0.06,
    rho_death = 0.7,
    phi_death = 2.2,
    alpha0   = 4, # beta-binomial intercept
    alpha1   = 0.8,    # beta-binomial slope
    kappa    = 2.2)

parameters =
  c(beta = true_pars[["R0"]] / popsize / (true_pars[["dur_early"]] + true_pars[["dur_progress"]]), # R0 = beta * P / mu
    gamma = 1 / true_pars[["dur_latent"]],
    nu_early = 1 / true_pars[["dur_early"]],
    mu_rec = (1 - true_pars[["ifr"]]) / true_pars[["dur_progress"]],
    mu_death = true_pars[["ifr"]] / true_pars[["dur_progress"]],
    rho_death = true_pars[["rho_death"]],
    phi_death = true_pars[["phi_death"]],
    alpha0 = true_pars[["alpha0"]],
    alpha1 = true_pars[["alpha1"]],
    kappa = true_pars[["kappa"]])

n_params <- length(parameters)
par_initializer = function() {
  priors$from_estimation_scale(priors$to_estimation_scale(parameters) + rnorm(n_params, 0, 0.1))
}


# Build Kernel ------------------------------------------------------------
# specify the kernel
mcmc_kern <-
  mcmc_kernel(
    parameter_blocks =
      list(parblock(
        pars_nat = c("beta", "gamma", "nu_early", "mu_rec", "mu_death", "rho_death", "phi_death", "alpha0", "alpha1", "kappa"),
        # par_est should have different names from pars_nat
        pars_est = c("R0_est", "dur_latent_est", "dur_early_est", "dur_progress_est", "ifr_est", "rho_death_est", "phi_death_est", "alpha0_est", "alpha1_est", "kappa_est"),
        priors = priors,
        alg = "mvnss",
        # alg = "mvnmh",
        sigma = diag(0.01, n_params),
        initializer = par_initializer,
        control =
          mvnss_control(stop_adaptation = 50000,
                                      scale_cooling = 0.85,
                                      scale_constant = 1,
                                      step_size = 0.25))),
    lna_ess_control = lna_control(bracket_update_iter = 5e3,
                                  joint_initdist_update = FALSE))


# Fit Model ---------------------------------------------------------------

# Old Fashion Way ---------------------------------------------------------


n_chains <- 4
thinning_interval <- 100
iterations <- 250000
multi_chain_stem_fit <- list()
multi_chain_stem_fit$n_chains <- n_chains
multi_chain_stem_fit$n_iterations <- (iterations - mcmc_kern$parameter_blocks[[1]]$control$stop_adaptation) / thinning_interval
# CHECK THIS
multi_chain_stem_fit$thinning_interval <- thinning_interval

multi_chain_stem_fit$stem_fit_list <- foreach(chain = 1:4,
                                              .packages = "stemr",
                                              .export = ls()) %dorng% {

                                                fit_stem(stem_object = stem_object,
                                                         method = "ode", # or "lna"
                                                         mcmc_kern = mcmc_kern,
                                                         iterations = iterations,
                                                         thinning_interval = thinning_interval,
                                                         print_progress = 1e3)
                                              }

write_rds(multi_chain_stem_fit,"mvnss_sim/mvnss_multi_chain_stem_fit.rds")
