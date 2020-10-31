library(tidyverse)


# Load Data ---------------------------------------------------------------
oc_data <- read_csv("/Users/damon/Documents/uci_covid_modeling/data/oc_data.csv",
                    col_types = cols(X1 = col_skip(), posted_date = col_date(format = "%Y-%m-%d"),
                                     new_cases = col_integer(),
                                     new_tests = col_integer(),
                                     new_deaths = col_integer(),
                                     new_deaths_calcat = col_skip())) %>%
  rename(date = posted_date)

lump_oc_data <- function(oc_data, time_interval_in_days, first_day = "0000-01-01", last_day = "9999-12-31") {
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

popsize <- 3068927L

dat <- lumped_oc_data %>%
  mutate(time = as.numeric((end_date - min(start_date) + 1L) / 7)) %>%
  select(time, cases = new_cases, deaths = new_deaths)

obs_times <- dat$time

oc_data %>%
  pivot_longer(-date) %>%
  group_by(name) %>%
  mutate(value = cumsum(value)) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  filter(date == "2020-03-14")

# Setup Stemr Model -------------------------------------------------------
strata <- NULL
compartments <- c("S", "E", "Ie", "Ip", "R", "D")

rates <-
  list(rate(rate = "beta_t * (Ie + 0.8 * Ip)", # individual level rate (unlumped)
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

init_infected <- sum(oc_data$new_cases[oc_data$date <= "2020-03-14"]) * 10
init_states <- c(popsize - init_infected, lengths(splitIndices(init_infected, 3)), 0, 0) %>%
  `names<-`(compartments)

state_initializer <-
  list(stem_initializer(
    init_states = init_states,
    fixed = F,
    prior = init_states / 10,
    dist = "dirmultinom"))



parameters <- c(gamma = 1, nu_early = 1, mu_rec = 0.94,
                mu_death = 0.06, rho_death = 0.7, phi_death = 2.2, alpha0 = 4,
                alpha1 = 0.8, kappa = 2.2)
constants <- c(t0 = 0)
tcovar <- data.frame(time = obs_times,
                     beta_t = (2.4 + 0.35 * sin(obs_times / 3.5) + rnorm(length(obs_times), 0, 0.1)) / 1e7,
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
    compile_lna = F,   # compile LNA functions
    messages = F       # don't print messages
  )


measurement_process <-
  stem_measure(emissions = emissions,
               dynamics = dynamics,
               data = dat)

stem_object <- make_stem(dynamics = dynamics, measurement_process = measurement_process)

stem_data <- simulate_stem(stem_object  = stem_object,
                           method       = "ode",
                           paths        = TRUE,
                           observations = T,
                           nsim         = 1,
                           census_times = c(0, obs_times))

stem_data
