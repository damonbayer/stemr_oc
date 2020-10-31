library(stemr)
library(tidybayes)
library(tidyverse)
source("~/Documents/stemr_oc/stemr_functions.R")
oc_prior <- read_rds("compare_beta_binomials/oc_prior.rds")
model_objects <-  read_rds("fixed_init_sim/model_objects.rds")
popsize <- model_objects$popsize

stan_params <- oc_prior %>%
  spread_draws(alpha_incid_0,
               alpha_incid_1,
               early_dur,
               IFR,
               latent_dur,
               prog_dur,
               prog_dur,
               sqrt_kappa_inv_incid,
               sqrt_phi_inv_death,
               log_R0,
               rho_death) %>%
  transmute(R0_est = log_R0,
            dur_latent_est = -log(latent_dur),
            dur_early_est = -log(early_dur),
            dur_progress_est = -log(prog_dur),
            ifr_est = logit(IFR),
            rho_death_est = logit(rho_death),
            phi_death_est = log(sqrt_phi_inv_death),
            alpha0_est = log(alpha_incid_0),
            alpha1_est = logit(alpha_incid_1),
            kappa_est = log(sqrt_kappa_inv_incid)) %>%
  as.matrix() %>%
  apply(1, from_estimation_scale) %>%
  t() %>%
  split(., row(.)) %>%
  map(~`names<-`(.,c("beta", "gamma", "nu_early", "mu_rec", "mu_death", "rho_death",
                     "phi_death", "alpha0", "alpha1", "kappa")))



# Prep Stemr --------------------------------------------------------------

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
                obstimes        = obs_times))
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

stemr_prior <- simulate_stem(stem_object = stem_object, method = "ode", nsim = 4000,
                             simulation_parameters = stan_params)


stemr_prior$
