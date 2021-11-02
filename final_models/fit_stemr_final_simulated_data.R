args <- commandArgs(trailingOnly=T)
my_seed <- as.numeric(args[1])

library(stemr)
library(extraDistr)
library(foreach)
library(doRNG)
library(doParallel)
library(tidyverse)
library(lubridate)
source("helper_functions.R")


popsize <- 3175692L
log_popsize <- log(popsize)


# Prepare OC Data ---------------------------------------------------------
first_day <- ymd("2020-03-30")
last_day <- ymd("2021-02-28")

if (Sys.info()[["sysname"]] == "Linux") {
  dat_raw <- read_csv("/data/homezvol2/bayerd/stemr_oc/final_models/data/simulated_data.csv")
} else if (Sys.info()[["sysname"]] == "Darwin") {
  dat_raw <- read_csv("final_models/data/simulated_data.csv")
}

dat <- dat_raw %>%
  select(-starts_with("sero"))

dat_seroprev <- dat_raw %>%
  select(time, ends_with("date"), starts_with("seroprev")) %>%
  drop_na()

obs_times <- dat$time
obs_times_seroprev <- dat_seroprev$time


# Build Model -------------------------------------------------------------
set.seed(my_seed)
strata <- NULL
compartments <- c("S", "I", "R", "D", "E")

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

if (args[1] == "double_initial_counts") {
  init_states <- c(S = 3168000, I = 5000, R = 1, D = 1, E = 2690)
  init_states_prior <- c(5.3252274084003, 0.619525040689332, -8.59062950948942,
                         -8.59044365315583, 0.12208870486078,
                         0.0641635068561201, 0.000237679033324254,
                         0.000237642407309249)
} else {
  init_states <- round(c(S = 3168000, I = 5000, R = 1, D = 1, E = 2690) * popsize / 3175692)
  init_states_prior <- c(6.02067492520774, 0.619153500651619, -7.89766815072691,
                         -7.89729647259588, 0.549603697818656,
                         0.0292277802326922, 9.8518767486955e-05,
                         9.84821737399636e-05)
}

state_initializer <-
  list(
    stem_initializer(
      init_states = init_states,
      fixed = F,
      prior = init_states_prior,
      dist = "sbln"
    )
  )

adjacency <- NULL

tcovar <- dat %>%
  select(time, tests, prop_deaths_reported) %>%
  left_join(dat_seroprev %>% select(time, seroprev_tests), by = "time") %>%
  as.matrix()


tmax <- max(obs_times, obs_times_seroprev)
constants <- c(t0 = 0)

foi_rw1 <- function(parameters, draws, log_pop = log_popsize) exp(c(log(parameters[["R0_init"]]), log(parameters[["R0_init"]]) + head(cumsum(draws), -1) * parameters[["sigma_R0"]]) - log_pop - log(parameters[["dur_infec"]]))

ifr_t_fn <- function(parameters, draws) c(parameters[["ifr_init"]], expit(logit(parameters[["ifr_init"]]) + head(cumsum(draws), -1) * parameters[["sigma_ifr"]]))

alpha_t_fn <- function(parameters, draws) c(parameters[["alpha_init"]], exp(log(parameters[["alpha_init"]]) + head(cumsum(draws), -1) * parameters[["sigma_alpha"]]))


tparam <-
  list(tpar(
    tparam_name = "beta_t",
    times = c(0, head(obs_times,-1)),
    n_draws = length(obs_times),
    draws2par = foi_rw1
  ),
  tpar(
    tparam_name = "ifr_t",
    times = c(0, head(obs_times, -1)),
    n_draws = length(obs_times),
    draws2par = ifr_t_fn
  ),
  tpar(
    tparam_name = "alpha_t",
    times = c(0, head(obs_times, -1)),
    n_draws = length(obs_times),
    draws2par = alpha_t_fn)
  )

parameters_initial <-
  c(
    R0_init = 1,
    gamma = 1,
    dur_infec = 1,
    rho_death = 0.7,
    phi_death = 2.2,
    alpha_init = 5,
    kappa = 2.2,
    ifr_init = 0.005,
    sigma_R0 = 0.08,
    sigma_ifr = 0.09,
    sigma_alpha = 0.05
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
          "kappa * (exp(alpha_t) * (E2I / popsize)) / (exp(alpha_t) * (E2I / popsize) + (1 - E2I / popsize))",
          "kappa * ((1 - E2I / popsize)) / ((E2I / popsize) + (1 - E2I / popsize))"
        ),
      incidence = TRUE,
      obstimes = obs_times
    ),
    emission(
      meas_var = "deaths",
      distribution = "negbinomial",
      emission_params = c("phi_death",
                          "rho_death * I2D * prop_deaths_reported"),
      incidence = TRUE,
      obstimes = obs_times
    ),
    emission(meas_var = "seroprev_cases",
             distribution = "binomial",
             emission_params = c("seroprev_tests", "R / (S + E + I + R)"),
             incidence = F,
             obstimes = obs_times_seroprev
    ))


measurement_process <- stem_measure(emissions = emissions,
                                    dynamics = dynamics,
                                    data = list(as.matrix(select(dat, time, cases)),
                                                as.matrix(select(dat, time, deaths)),
                                                as.matrix(select(dat_seroprev, time, seroprev_cases))))

stem_object <-
  make_stem(dynamics = dynamics, measurement_process = measurement_process)


# Set Up Parameters -------------------------------------------------------
to_estimation_scale = function(params_nat) {
  c(R0_init_est = log(params_nat[["R0_init"]]), # log(R0_init)
    dur_latent_est = log(params_nat[["gamma"]]), # -log(dur_latent)
    dur_infec_est = log(params_nat[["dur_infec"]]), # log(dur_infec)
    rho_death_est = logit(params_nat[["rho_death"]]), # logit(rho_death)
    phi_death_est = -0.5 * log(params_nat[["phi_death"]]), # -0.5 * log(phi_death)
    alpha_init_est = log(params_nat[["alpha_init"]]), # log(alpha)
    kappa_est = -0.5 * log(params_nat[["kappa"]]), # -0.5 * log(kappa)
    ifr_init_est = logit(params_nat[["ifr_init"]]),
    sigma_R0_est = log(params_nat[["sigma_R0"]]),
    sigma_ifr_est = log(params_nat[["sigma_ifr"]]),
    sigma_alpha_est = log(params_nat[["sigma_alpha"]]))
}

# Need to account for sigma
from_estimation_scale = function(params_est) {
  c(R0_init = exp(params_est[["R0_init_est"]]),
    gamma = exp(params_est[["dur_latent_est"]]),
    dur_infec = exp(params_est[["dur_infec_est"]]),
    rho_death = expit(params_est[["rho_death_est"]]),
    phi_death = exp(-2 * params_est[["phi_death_est"]]),
    alpha_init = exp(params_est[["alpha_init_est"]]),
    kappa = exp(-2 * params_est[["kappa_est"]]),
    ifr_init = expit(params_est[["ifr_init_est"]]),
    sigma_R0 = exp(params_est[["sigma_R0_est"]]),
    sigma_ifr = exp(params_est[["sigma_ifr_est"]]),
    sigma_alpha = exp(params_est[["sigma_alpha_est"]]))
}


if (args[1] == "lower_R0") {
  R0_init_est_mean <- log(0.5)
  R0_init_est_sd <- 0.4
} else {
  R0_init_est_mean <- 0
  R0_init_est_sd <- 0.225
}
dur_latent_est_mean <- log(3/7)
dur_latent_est_sd <- 0.25
if (args[1] == "longer_infectious_period") {
  dur_infec_est_mean <- log(5/7)
  dur_infec_est_sd <- 0.2
} else {
  dur_infec_est_mean <- log(7/7)
  dur_infec_est_sd <- 0.225
}
rho_death_est_shape1 <- 200
rho_death_est_shape2 <- 20
phi_death_est_rate <- 1
if (args[1] == "lower_alpha") {
  alpha_init_est_mean <- 2
  alpha_init_est_sd <- 0.5
} else {
  alpha_init_est_mean <- 4
  alpha_init_est_sd <- 0.5
}
kappa_est_rate <- 1
if (args[1] == "higher_ifr") {
  ifr_init_est_shape1 <- 60
  ifr_init_est_shape2 <- 9000
} else {
  ifr_init_est_shape1 <- 30
  ifr_init_est_shape2 <- 10000
}
sigma_R0_est_mean <- 0.05
sigma_R0_est_sd <- 0.01
sigma_ifr_est_mean <- 0.08
sigma_ifr_est_sd <- 0.01
sigma_alpha_est_mean <- 0.05
sigma_alpha_est_sd <- 0.01

logprior =
  function(params_est) {
    sum(
      dnorm(params_est["R0_init_est"], mean = R0_init_est_mean, sd = R0_init_est_sd, log = TRUE), # log(R0)
      dnorm(-params_est["dur_latent_est"],  mean = dur_latent_est_mean, sd = dur_latent_est_sd, log = TRUE), # -log(dur_latent)
      dnorm(params_est["dur_infec_est"], mean = dur_latent_est_mean, sd = dur_infec_est_sd, log = T),
      dbeta(expit(params_est["rho_death_est"]), shape1 = rho_death_est_shape1, shape2 = rho_death_est_shape2, log = TRUE) + params_est["rho_death_est"] - 2 * log(exp(params_est["rho_death_est"]) + 1) , # logit(rho_death)
      dexp(exp(params_est["phi_death_est"]), rate = phi_death_est_rate, log = TRUE) +  params_est["phi_death_est"], # -0.5 * log(phi_death)
      dtnorm(exp(params_est["alpha_init_est"]), mean = alpha_init_est_mean, sd = alpha_init_est_sd, a = 0, log = TRUE) + params_est["alpha_init_est"], # log(alpha)
      dexp(exp(params_est["kappa_est"]), rate = kappa_est_rate, log = T) +  params_est["kappa_est"], # -0.5 * log(kappa)
      dbeta(expit(params_est["ifr_init_est"]), shape1 = ifr_init_est_shape1, shape2 = ifr_init_est_shape2, log = TRUE) + params_est["ifr_init_est"] - 2 * log(exp(params_est["ifr_init_est"]) + 1), # logit(ifr_init)
      dtnorm(exp(params_est["sigma_R0_est"]), mean = sigma_R0_est_mean, sd = sigma_R0_est_sd, a = 0, log = TRUE) + params_est["sigma_R0_est"], # log(sigma_R0)
      dtnorm(exp(params_est["sigma_ifr_est"]), mean = sigma_ifr_est_mean, sd = sigma_ifr_est_sd, a = 0, log = TRUE) + params_est["sigma_ifr_est"], # log(sigma_ifr)
      dtnorm(exp(params_est["sigma_alpha_est"]), mean = sigma_alpha_est_mean, sd = sigma_alpha_est_sd, a = 0, log = TRUE) + params_est["sigma_alpha_est"]  # log(sigma_alpha)
    )
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

initial_values <- structure(c(
  1.15122629457544, 1.25649794322914, 1.25854844884716,
  1.13625227396481, 1.19172451895629, 1.13358797057164, 1.20210650934144,
  1.21104178098218, 1.21785410365415, 1.17996391707491, 3.26886100565582,
  3.28873423357374, 4.19984957896454, 3.63920069282901, 2.58812006226809,
  2.85560586842, 2.93918938772882, 2.89413498597635, 3.39743218804016,
  3.29489290668073, 0.733227364428071, 0.878995982668305, 0.915161896050096,
  0.693661919610873, 0.773411520598913, 0.744090364904359, 0.668547886505264,
  0.729748672093996, 0.805771174011239, 0.866834917362858, 0.908879913114176,
  0.881619040323954, 0.904765436536343, 0.921125924973078, 0.900928615255703,
  0.904508437476229, 0.919346618128245, 0.949085287946538, 0.932171372511707,
  0.920460459654317, 239.81835678338, 194.711969107417, 321.234944256983,
  167.390327769356, 268.365686216967, 290.331126862661, 379.534891337098,
  329.162951219892, 370.817715894835, 430.527655839715, 3.97448772455812,
  3.88094267180023, 4.15444637363217, 4.06618721684312, 4.3652458474725,
  4.48262593448254, 4.12357369870084, 4.18679285376786, 4.25545452495909,
  4.353397657496, 963873286.749679, 908171152.5821, 463789024.064355,
  571569180.931776, 3394038986.96009, 9230944673.85418, 758890002.058136,
  26080771.9901199, 183432850.789571, 184512343.099759, 0.00284695022744819,
  0.00296267414887485, 0.00275400163393943, 0.00269296584141162,
  0.00313301384139115, 0.00219533670931258, 0.0030813177460639,
  0.00320714979836539, 0.00309473909298774, 0.00319770922680797,
  0.0644831052493849, 0.0620133649270813, 0.0596902410381136, 0.0612011550947626,
  0.0740378586123304, 0.0656656217535556, 0.060420525489056, 0.0704745980204676,
  0.0686270922370653, 0.061976985547675, 0.07688094836616, 0.0827196051282812,
  0.0999611959929428, 0.0810838335154676, 0.0758894473267132, 0.0995305474924498,
  0.0838463076508917, 0.0712754804205866, 0.0924524871785279, 0.0771542855244768,
  0.0773990152801426, 0.0798115737254148, 0.0668665462565269, 0.0759076736314645,
  0.0687903044167596, 0.0803147651214962, 0.0712809933606284, 0.0723179280827186,
  0.0672260859059286, 0.0588415162180674
),
.Dim = 10:11,
.Dimnames = list(
  NULL, c(
    "R0_init", "gamma", "dur_infec", "rho_death", "phi_death",
    "alpha_init", "kappa", "ifr_init", "sigma_R0", "sigma_ifr", "sigma_alpha"
  )
)
)[my_seed, ]

custom_par_initializer <- function() initial_values

to_human_scale <- function(params_est) {
  c(R0_init = exp(params_est[["R0_init_est"]]),
    dur_latent = exp(-params_est[["dur_latent_est"]]),
    dur_infec = exp(params_est[["dur_infec_est"]]),
    rho_death = stemr::expit(params_est[["rho_death_est"]]),
    phi_death = exp(-2 * params_est[["phi_death_est"]]),
    alpha = exp(params_est[["alpha_init_est"]]),
    kappa = exp(-2 * params_est[["kappa_est"]]),
    ifr_init = stemr::expit(params_est[["ifr_init_est"]]),
    sigma_R0 = exp(params_est[["sigma_R0_est"]]),
    sigma_ifr = exp(params_est[["sigma_ifr_est"]]),
    sigma_alpha = exp(params_est[["sigma_alpha_est"]]))
}

mcmc_kern <-
  mcmc_kernel(
    parameter_blocks =
      list(parblock(
        pars_nat = c("R0_init", "gamma", "dur_infec", "rho_death", "phi_death", "alpha_init", "kappa", "ifr_init", "sigma_R0", "sigma_ifr", "sigma_alpha"),
        pars_est = c("R0_init_est", "dur_latent_est", "dur_infec_est", "rho_death_est", "phi_death_est", "alpha_init_est", "kappa_est", "ifr_init_est", "sigma_R0_est", "sigma_ifr_est", "sigma_alpha_est"),
        priors = priors,
        alg = "mvnss",
        sigma = diag(0.01, n_params),
        initializer = custom_par_initializer,
        control =
          mvnss_control(stop_adaptation = 50000,
                        scale_cooling = 0.85,
                        scale_constant = 1,
                        step_size = 0.25))),
    lna_ess_control = lna_control(bracket_update_iter = 5e3,
                                  joint_initdist_update = FALSE),
    tparam_ess_control = tpar_control(bracket_update_iter = 1e4))

n_chains <- 1
thinning_interval <- 3500
iterations <- 12500000

res <- fit_stem(stem_object = stem_object,
                method = "ode", # or "lna"
                mcmc_kern = mcmc_kern,
                iterations = iterations,
                thinning_interval = thinning_interval,
                print_progress = 1e3)

multi_chain_stem_fit <- list()
multi_chain_stem_fit$n_iterations <- iterations
multi_chain_stem_fit$n_chains <- n_chains
multi_chain_stem_fit$thinning_interval <- thinning_interval
multi_chain_stem_fit$stem_fit_list <- res
multi_chain_stem_fit$to_human_scale <- to_human_scale
multi_chain_stem_fit$data <- left_join(dat, select(dat_seroprev, -ends_with("date")))
write_rds(multi_chain_stem_fit, paste0("res_fit_stemr_final_models_simulated_data_", my_seed, "_", format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), ".rds"))
