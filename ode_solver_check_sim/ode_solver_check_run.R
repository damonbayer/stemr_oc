source('stemr_functions.R')
source("~/Documents/uci_covid_modeling/code/SEIeIpRD/plot_functions.R")
library(tidyverse)
library(tidybayes)
library(fs)
library(rstan)
library(stemr)
theme_set(cowplot::theme_cowplot())
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


# Setup Stan --------------------------------------------------------------
model_objects <- read_rds("fixed_init_sim/model_objects.rds")
model_objects_priors_only <- model_objects
model_objects_priors_only$priors_only <- T
model_objects_priors_only$forecast_in_days <- 0

control_list <- list(adapt_delta = 0.99,
                     max_treedepth = 20)


# Setup Stemr -------------------------------------------------------------
# Read Data ---------------------------------------------------------------
dat <- read_rds("~/Documents/uci_covid_modeling/code/SEIeIpRD/oc/2020-08-12_2020-09-16/model_objects.rds")$lumped_ochca_covid

init_states <- c(S = 3012564L, E = 14866L, Ie = 20907L, Ip = 20590L, R = 0L, D = 0L)
popsize <- sum(init_states)

if (popsize != model_objects$popsize) {stop("popsizes disagree")}

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


# Run Stan ----------------------------------------------------------------
stan_results <- stan(file = "~/Documents/stemr_oc/fixed_init_sim/SEIeIpRD.stan",
                 data = model_objects_priors_only,
                 algorithm = "Fixed_param", iter = 1,
                 seed = 0,
                 chains = 1,
                 control = control_list)

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

stan_params <- stan_results %>%
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
  pivot_longer(everything()) %>%
  deframe() %>%
  from_estimation_scale()

stemr_results <- simulate_stem(stem_object = stem_object, method = "ode", simulation_parameters = list(stan_params))



# Compare -----------------------------------------------------------------
stemr_epi_curves <- stemr_results$natural_paths[[1]] %>%
  as_tibble() %>%
  pivot_longer(-time) %>%
  mutate(model = "stemr")

init_fracs <- spread_draws(stan_results, init_fracs[i])$init_fracs

stan_epi_curves <- prep_epi_curves(stan_results, model_objects) %>%
  pivot_wider( values_from = epi_curves, names_from = g_text) %>%
  select(.draw, t, S, E, IE, IP, D) %>%
  mutate(R = 1 - sum(S, E, IE, IP, D)) %>%
  pivot_longer(-c(.draw, t)) %>%
  rename(time = t) %>%
  ungroup() %>%
  mutate(time = as.numeric(time - min(time) + 3) / 7,
         name = str_to_title(name)) %>%
  filter(time %in% unique(stemr_epi_curves$time)) %>%
  filter(name %in% unique(stemr_epi_curves$name)) %>%
  mutate(value = value * model_objects$popsize) %>%
  select(time, name, value) %>%
  mutate(name = fct_inorder(name),
         model = "stan") %>%
  bind_rows(tibble(S = 1 - init_fracs[1],
                   E = exp(log(init_fracs[1]) + log(1 - init_fracs[2])),
                   Ie = exp(log(init_fracs[1]) + log(init_fracs[2]) + log(init_fracs[3])),
                   Ip = exp(log(init_fracs[1]) + log(init_fracs[2]) + log(1-init_fracs[3])),
                   R = 0,
                   D = 0) %>%
              pivot_longer(everything()) %>%
              mutate(value = value * model_objects$popsize) %>%
              mutate(time = 0, model = "stan"), .)

combined_epi_curves <- bind_rows(stemr_epi_curves, stan_epi_curves)

write_rds(combined_epi_curves, "ode_solver_check_sim/combined_epi_curves.rds")

combined_epi_curves %>%
  ggplot(aes(time, value, group = model, color = model)) +
  geom_line(size = 1.25, alpha = 0.5) +
  facet_wrap(. ~ name, scales = "free_y")
