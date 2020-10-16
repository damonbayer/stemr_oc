library(tidyverse)
library(stemr)
library(extraDistr)
theme_set(cowplot::theme_cowplot())

# Read Data ---------------------------------------------------------------
dat <- read_rds("~/Documents/uci_covid_modeling/code/SEIeIpRD/oc/2020-08-12_2020-09-16/model_objects.rds")$lumped_ochca_covid

model_objects_path <- NULL
if (is.null(model_objects_path)) stop("please download this file and update model_objects_path\nhttps://github.com/vnminin/uci_covid_modeling/raw/master/code/SEIeIpRD/oc/2020-08-12_2020-09-16/model_objects.rds")


dat <- read_rds(model_objects_path)
# popsize <- 3175692L
# init_states <- c(S = 3012564L, E = 14866L, Ie = 20907L, Ip = 20590L, R = 106242L, D = 522L)

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

# parameters <- numeric(10); names(parameters) <- c("beta", "gamma", "nu_early", "mu_rec", "mu_death", "rho_death", "phi_death", "alpha0", "alpha1", "kappa")
parameters <- c(beta = 1.40951517189436e-07, gamma = 1.00205484766246, nu_early = 1.00047239762145,
                mu_rec = 0.991220435483575, mu_death = 0.00590274654390186, rho_death = 0.822026217472393,
                phi_death = 2.08767492297483, alpha0 = 57.1282936397807, alpha1 = 0.79420146357785,
                kappa = 2.08114643628088)
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

tmp <- simulate_stem(stem_object = stem_object, method = "ode")
