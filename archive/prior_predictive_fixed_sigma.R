library(stemr)
library(extraDistr)
library(foreach)
library(doRNG)
library(doParallel)
library(tidyverse)
library(tidybayes)
source('~/Documents/stemr_oc/stemr_functions.R')
my_sigma <- 0.2

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

lumped_oc_data <- lump_oc_data(oc_data, time_interval_in_days = 7, first_day = "2020-03-14", last_day = "2020-10-16")

dat <- lumped_oc_data %>%
  mutate(time = as.numeric((end_date - min(start_date) + 1L) / 7)) %>%
  select(time, cases = new_cases, deaths = new_deaths, tests = new_tests)

obs_times <- dat$time

popsize <- 3068927L
log_popsize <- log(popsize)


# Build Model -------------------------------------------------------------
# set.seed(12511)
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
  sum(oc_data$new_cases[oc_data$date <= "2020-03-14"]) * 30

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
    R0_init = 2,
    gamma = 1,
    nu_early = 1,
    mu_rec = 0.94,
    mu_death = 0.06,
    rho_death = 0.7,
    phi_death = 2.2,
    alpha0 = 4,
    alpha1 = 0.8,
    kappa = 2.2
  )

constants <- c(t0 = 0, sigma = my_sigma)
tmax <- max(obs_times)

foi_rw1 <- function(parameters, draws, log_pop = log_popsize) {
  log_R0_t <- numeric(length = length(draws))
  log_R0_t[1] <- log(parameters["R0_init"])

  for (t in 2:(length(log_R0_t))) {
    log_R0_t[t] <- log_R0_t[t - 1] + draws[t - 1] * constants[["sigma"]]
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

from_estimation_scale = function(params_est) {
  c(R0_init = exp(params_est[["R0_init_est"]]),
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

set.seed(200)
simulation_parameters_est <- rerun(.n = 2000, {
  c(R0_init_est = rnorm(1, -0.2554128198465173693599, 0.7),
    dur_latent_est = -rnorm(1, 0, 0.22),
    dur_early_est = -rnorm(1, 0, 0.22),
    dur_progress_est = -rnorm(1, 0, 0.22),
    ifr_est = logit(rbeta(1, 1.5, 200)),
    rho_death_est = logit(rbeta(1, 8, 2)),
    phi_death_est = -2 * log(rexp(1, 1)),
    alpha0_est = log(rtnorm(1, mean = 4, sd = 2, a = 0)),
    alpha1_est = logit(rbeta(1, 3, 1)),
    kappa_est = -2 * log(rexp(1, 1)^2))
})

tparam_draws <- map(simulation_parameters_est, ~list(foi_rw1(parameters = from_estimation_scale(.), draws = rnorm(31))))
# simulation_parameters <- map(simulation_parameters_est, ~c(from_estimation_scale(.), tail(stem_object$dynamics$parameters, 6)))
simulation_parameters <- map(simulation_parameters_est, from_estimation_scale)

prior_pred <- simulate_stem(stem_object = stem_object,
                            nsim = 2000,
                            simulation_parameters = simulation_parameters,
                            tparam_draws = tparam_draws,
                            method = "ode", paths = F)


{
  beta_t <- map(tparam_draws, pluck(1)) %>% unlist() %>% array(dim = c(31, 1, length(tparam_draws)))
  nu_early <- map_dbl(simulation_parameters, "nu_early")
  mu_death <- map_dbl(simulation_parameters, "mu_death")
  mu_rec <- map_dbl(simulation_parameters, "mu_rec")
}

rt_tmp <-
  exp(log(beta_t) + log_popsize + rep(log(1 / nu_early + 0.8 / (mu_rec + mu_death)), each = dim(beta_t)[2])) %>%
  gather_array(value = value, time, var, .iteration) %>%
  as_tibble() %>%
  select(-var) %>%
  mutate(time = time - 1) %>%
  mutate(.chain = 1) %>%
  mutate(.draw = tidybayes:::draw_from_chain_and_iteration_(.chain, .iteration)) %>%
  left_join(prior_pred$natural_paths %>%
              do.call(what = rbind, args = .) %>%
              as_tibble() %>%
              mutate(.chain = 1,
                     .iteration = rep(1:2000, each = length(unique(time))),
                     .draw = tidybayes:::draw_from_chain_and_iteration_(.chain, .iteration)) %>%
              mutate(prop_s = S / popsize) %>%
              select(time, .iteration, .chain, .draw, prop_s)) %>%
  mutate(Rt = value * prop_s)

# Plot betas
# rt_tmp %>%
#   select(time, value) %>%
#   group_by(time) %>%
#   median_qi(.width = c(0.5, 0.8, 0.95)) %>%
#   ggplot(aes(time, value, ymin = .lower, ymax = .upper)) +
#   geom_lineribbon() +
#   geom_hline(yintercept = 0.2, color = "red") +
#   scale_fill_brewer() +
#   theme_minimal()

rt_plot_020 <- {
ggplot() +
  geom_lineribbon(
    data = rt_tmp %>%
      select(time, Rt) %>%
      group_by(time) %>%
      median_qi(.width = c(0.5, 0.8, 0.95)),
    mapping = aes(time, Rt, ymin = .lower, ymax = .upper)) +
  geom_line(data = rt_tmp %>%
              filter(.iteration %in% sample(2000, 12)),
            mapping = aes(time, Rt, group = .iteration)) +
  geom_hline(yintercept = c(0.2, 4), color = "red") +
  scale_fill_brewer() +
  theme_minimal() +
  scale_y_log10() +
  theme(legend.position = "none") +
  ggtitle(str_c("Sigma = ", my_sigma))
  }

rt_plot_005 + scale_y_log10(limits = c(0.05, 5))

rt_plot_020 + scale_y_log10(limits = c(0.025, 40))

cowplot::plot_grid(rt_plot_005 + scale_y_log10(limits = c(0.025, 40)),
                   rt_plot_010 + scale_y_log10(limits = c(0.025, 40)),
                   rt_plot_015 + scale_y_log10(limits = c(0.025, 40)),
                   rt_plot_020 + scale_y_log10(limits = c(0.025, 40)), align = "hv")
