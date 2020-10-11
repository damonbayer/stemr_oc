#' ---
#' title: "stemr: Baysian Inference for Stochastic Epidemic Models via the Linear Noise Approximation"
#' subtitle: Simulating from and fitting an SIR model
#' author: "Jonathan Fintzi, Jon Wakefield, and Vladimir N. Minin"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{stemr: Baysian Inference for Stochastic Epidemic Models via the Linear Noise Approximation}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' ---
#'
## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  message = FALSE,
  comment = "#>"
)

#'
#' # Overview
#'
#' This vignette demonstrates the basic functionalities of the  `stemr` package.
#' Broadly, the package simulates and fits stochastic epidemic models to partially
#' observed incidence and prevalence data, and implements the Bayesian data
#' augmentation framework presented in Fintzi, Wakefield, and Minin (2019).
#' Simulation can be conducted using ordinary  differential equations (ODEs),
#' Gillespie's direct algorithm (Gillespie, 1976), and using a restarting version
#' of the linear noise approximation (LNA; Fintzi et al., 2019). Inference is
#' conducted using ODEs or the LNA. This vignette, in particular, demonstrates how
#' to simulating partially observed incidence data from an outbreak with SIR
#' dyanmics and fit an SIR model to the data via the linear noise approximation.
#' This corresponds to the procedure used in the coverage simulation in Fintzi et
#' al., (2019).
#'
#' # Installing and loading the `stemr` package
#'
#' To install the `stemr` package, clone this repository and build the package from sources. There are two important things to note. First, the package depends on having version 0.4.1 of the `Ryacas` package installed. This can be accomplished, using the `devtools` package: `devtools::install_version("Ryacas", "0.4.1")`. It is also critical that the `stemr` package is installed without byte compilation. See [this page](https://support.rstudio.com/hc/en-us/articles/200486518-Customizing-Package-Build-Options) for how to do this. You should be able to rebuild in the usual way once you clone the package repo and install the other dependencies (odeintr, MASS, extraDistr, stats, ggplot2, cowplot, Rcpp, RcppArmadillo, and BH).
#'
## ---- include=FALSE------------------------------------------------------
require(ggplot2)
require(cowplot)
set.seed(12511)

#'
#' # Basic example: partially observed incidence from an outbreak with SIR dynamics
#'
#' As a basic example, we will simulate an outbreak with SIR dynamics and generate
#' negative binomial incidence counts with a mean case detection rate of 0.5. This
#' corresponds to the coverage simulation conducted in Fintzi et al. (2019).
#' Briefly, the SIR model describes the time-evolution of an outbreak homogeneously
#' mixing population where each individual exists in one of three states -
#' susceptible (S), infected (I), and recovered (R). Infection implies that an
#' individual is infectious with no latent period, while recovery confers lifelong
#' immunity. The waiting times between infection and recovery events are taken to
#' be exponentially distributed and the rates of state transition change every time
#' the population jumps to a different state. Hence, the transmission process is a
#' Markov jump process. The rates at which infection and recovery events take place
#' are $$\lambda_{SI} = \beta S I,\ \text{and}\ \lambda_{IR} = \mu I,$$
#' where $\beta$ is the per-contact rate of infection, $\mu$ is the recovery rate,
#' and $S$ and $I$ denote the numbers of susceptible and infectious individuals.
#' The basic reproduction number under SIR dynamics is $R_0 = \beta P/\mu$, where
#' $P=S+I+R$ is the population size.
#'
#' We initialize parameters and instatiate the SIR model in the code block below.
#' The functions in the `stemr` package are documented and can be accessed in the usual way,
#' e.g., `help(rate)`.
#'
## ---- echo = TRUE, warnings = FALSE--------------------------------------
library(stemr)
popsize = 1e4 # population size
# popsize <- 3175692L

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
        kappa    = 2.2)     # beta-binomial overdispersion

# initialize model compartments and rates
strata <- NULL # no strata
compartments <- c("S", "E", "Ie", "Ip", "R", "D")

# set the parameter values - must be a named vector
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

# rates initialized as a list of rate lists
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

# list used for simulation/inference for the initial state, initial counts fixed.
# state initializer a list of stem_initializer lists.
state_initializer <-
  list(stem_initializer(
          init_states = c(S = popsize-10, E = 10, Ie = 0, Ip = 0, R = 0, D = 0), # must match compartment names
          fixed = T)) # initial state fixed for simulation, we'll change this later

# declare the initial time to be constant
constants <- c(t0 = 0, tests = 1000)
tmax <- 40
tests <- 1000
# tcovar <- data.frame(time = 1:tmax,
                   # tests = rep(seq(1000, 4000, by = 1000),each = 10))


# compile the model
dynamics <-
      stem_dynamics(
            rates = rates,
            tmax = tmax,
            parameters = parameters,
            state_initializer = state_initializer,
            # tcovar = tcovar,
            compartments = compartments,
            constants = constants,
            compile_ode = T,   # compile ODE functions
            compile_rates = T, # compile MJP functions for Gillespie simulation
            compile_lna = T,   # compile LNA functions
            messages = F       # don't print messages
      )

# list of emission distribution lists (analogous to rate specification)
emissions <-
  list(emission(meas_var = "cases", # transition or compartment being measured (S->I transitions)
                distribution    = "betabinomial",        # emission distribution
                emission_params =
                  c("tests",
                    "kappa * (alpha0 * (Ie2Ip / popsize) ^ alpha1) / (alpha0 * (Ie2Ip / popsize) ^ alpha1 + (1 - Ie2Ip/popsize) ^ alpha1)",
                    "kappa * ((1 - Ie2Ip/popsize) ^ alpha1) / (alpha0 * (Ie2Ip / popsize) ^ alpha1 + (1 - Ie2Ip/popsize) ^ alpha1)"), # distribution pars, here overdispersion and mean
                incidence       = TRUE,                  # is the data incidence
                obstimes        = seq(1, tmax, by =1)), # vector of observation times
       emission(meas_var = "deaths",
                distribution = "negbinomial",
                emission_params = c("phi_death",
                                    "rho_death * Ip2D"),
                incidence = T,
                obstimes        = seq(1, tmax, by =1))) # vector of observation times)

# compile the measurement process
measurement_process <-
  stem_measure(emissions = emissions,
               dynamics  = dynamics,
               messages  = F)

# put it all together into a stochastic epidemic model object
stem_object <-
  make_stem(dynamics = dynamics,
       measurement_process = measurement_process)

#'
#' ## Simulating an outbreak and data
#'
#' Having compiled the model, we can simulate an outbreak and incidence (or
#' prevalence) data using Gillespie's direct algorith to simulate a MJP path, with
#' the LNA, or deterministically with ODEs. Since we specified that the data should
#' consist of incidence counts in instatiating the model, the simulation function
#' will produce incidence data. The prevalence and incidence curves are plotted
#' below.
#'
## ----sim_paths, echo = TRUE----------------------------------------------
sim_mjp <- simulate_stem(stem_object = stem_object, method = "gillespie", full_paths = T)
sim_lna <- simulate_stem(stem_object = stem_object, method = "lna", lna_method = "approx")
sim_ode <- simulate_stem(stem_object = stem_object, method = "ode")

#'
## ----plot_sims, echo = FALSE, fig.width=8, fig.height=4------------------
sim_paths =
    expand.grid(time = 0:tmax,
                Method = c("Gillespie", "LNA", "ODE"),
                Compartment = c("S","I","R","S2I","I2R"),
                Type = c("Prevalence","Incidence"))
sim_paths =
  sim_paths[!((
    sim_paths$Compartment %in% c("S", "I", "R") &
      sim_paths$Type == "Incidence"
  ) |
    (
      sim_paths$Compartment %in% c("S2I", "I2R") &
        sim_paths$Type == "Prevalence"
    )
  ), ]
sim_paths$Compartment = factor(sim_paths$Compartment, levels = c("S", "I", "R", "S2I", "I2R"))
sim_paths = sim_paths[with(sim_paths, order(Method, Compartment, Type, time)), ]
sim_paths$Count =
  c(
    sim_mjp$paths[[1]][, -1],
    sim_lna$natural_paths[[1]][, -1],
    sim_lna$paths[[1]][, -1],
    sim_ode$natural_paths[[1]][, -1],
    sim_ode$paths[[1]][, -1]
  )

mjp_prev =
  data.frame(time = sim_mjp$full_paths[[1]][,1],
             Compartment = rep(c("S","I","R"), each = nrow(sim_mjp$full_paths[[1]])),
             Count = c(sim_mjp$full_paths[[1]][,3:5]))

mjp_counts =
  ggplot(mjp_prev, aes(x = time, y = Count,
                       colour = Compartment,
                       group = Compartment)) +
  geom_step() +
  theme_minimal() +
  scale_color_brewer(type = "qual", palette = 6) +
  scale_y_continuous(trans = "sqrt",
                     breaks = c(0,50, 250, 1000, 2.5e3, 5e3,7.5e3,1e4),
                     expand = c(0,0)) +
  labs(title = "Compartment counts", subtitle = "MJP")

mjp_incid =
  ggplot(subset(sim_paths, Method == "Gillespie" & Type == "Incidence"),
         aes(x = time, y = Count,
             colour = Compartment,
             group = Compartment)) +
  geom_point() +
  theme_minimal() +
  scale_color_brewer("Transition", type = "qual", palette = 2) +
  labs(title = "Incident transition events",subtitle = "MJP")

lna_prev =
  ggplot(subset(sim_paths, Method == "LNA" & Type == "Prevalence"),
         aes(x = time, y = Count,
             colour = Compartment,
             group = Compartment)) +
  geom_line(linetype = 2) +
  theme_minimal() +
  scale_color_brewer("Compartment", type = "qual", palette = 6) +
  scale_y_continuous(trans = "sqrt",
                     breaks = c(0,50, 250, 1000, 2.5e3, 5e3,7.5e3,1e4),
                     expand = c(0,0)) +
  labs(title = "", subtitle = "LNA")

lna_incid =
  ggplot(subset(sim_paths, Method == "LNA" & Type == "Incidence"),
         aes(x = time, y = Count,
             colour = Compartment,
             group = Compartment)) +
  geom_point(shape = 2) +
  theme_minimal() +
  scale_color_brewer("Transition", type = "qual", palette = 2) +
  labs(title = "",subtitle = "LNA")

ode_prev =
  ggplot(subset(sim_paths, Method == "ODE" & Type == "Prevalence"),
         aes(x = time, y = Count,
             colour = Compartment,
             group = Compartment)) +
  geom_line(linetype = 3) +
  theme_minimal() +
  scale_color_brewer("Compartment", type = "qual", palette = 6) +
  scale_y_continuous(trans = "sqrt",
                     breaks = c(0,50, 250, 1000, 2.5e3, 5e3,7.5e3,1e4),
                     expand = c(0,0)) +
  labs(title = "", subtitle = "ODE")

ode_incid =
  ggplot(subset(sim_paths, Method == "ODE" & Type == "Incidence"),
         aes(x = time, y = Count,
             colour = Compartment,
             group = Compartment)) +
  geom_point(shape = 3) +
  theme_minimal() +
  scale_color_brewer("Transition", type = "qual", palette = 2) +
  labs(title = "", subtitle = "ODE")

cowplot::plot_grid(
  mjp_counts, lna_prev, ode_prev, mjp_incid, lna_incid, ode_incid
)

#'
#' ## Fitting the model to data
#'
#' We'll keep the dataset simulated under the MJP (shown below) and fit an SIR
#' model to the data.
#'
## ---- plot_dat-----------------------------------------------------------
ggplot(data = as.data.frame(sim_mjp$datasets[[1]]),
       aes(x=time, y = cases)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Week", y = "Count", title = "Observed Incidence")

#'
#' We'll need to recompile the measurement process with the simulated data fed to
#' the data argument of the `stem_measure` function since the data was not present
#' when the measurement process was originally compiled. There is no need to recompile
#' the model dynamics object if nothing about the model dynamics changes. In this
#' case however, we'll also do inference on the initial compartment volumes, so
#' we'll recompile the dynamics with a new state initializer.
#'
## ----recompile_meas, echo = TRUE-----------------------------------------

state_initializer <-
  list(stem_initializer(
    init_states = c(S = popsize-10, E = 10, Ie = 0, Ip = 0, R = 0, D = 0), # must match compartment names
    fixed = F,
    prior = c(popsize, 10, 0, 0, 0, 0) / 10
    )) # initial state fixed for simulation, we'll change this later

dynamics <-
  stem_dynamics(
    rates = rates,
    tmax = tmax,
    parameters = parameters,
    state_initializer = state_initializer,
    compartments = compartments,
    constants = constants,
    # tcovar = tcovar,
    compile_ode = T,   # compile ODE functions
    compile_rates = F, # compile MJP functions for Gillespie simulation
    compile_lna = T,   # compile LNA functions
    messages = F       # don't print messages
  )

measurement_process <-
  stem_measure(emissions = emissions,
               dynamics = dynamics,
               data = sim_mjp$datasets[[1]])

stem_object <- make_stem(dynamics = dynamics, measurement_process = measurement_process)

#'
#' In order to perform inference, we'll need to specify a function for
#' transforming the model parameters from their natural scale to the estimation
#' scale on which the MCMC explores the posterior, a function for transforming
#' parameters on their estimation scale to the natural scale on which they enter
#' the model dynamics and measurement process, and a function that returns the
#' log prior. We'll parameterize the MCMC estimation scale in terms of the log
#' basic reproduction number, log infectious period duration rate, and
#' unconstrained beta-binomial hyperparameters. The functions are specified as
#' follows and placed into a list of functions (note that it is critical that
#' the function signatures follow the specification given below):
#'

## ----priors_est_scale_fcns, echo = TRUE----------------------------------

### Parameterization in terms of log(R0) and log(mu)
## Priors for log(R0), log(mu), logit(rho), phi
# Parameters (natural scale): beta, mu, rho, phi
# Parameters (estimation scale): log(beta * N / mu), log(mu), logit(rho), log(phi)

# function to take params_nat and return params_est
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
#DONE
# function to take params_est and return params_nat
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
#DONE
# calculate the log prior density. note the jacobian for phi
logprior =
      function(params_est) {
        sum(dnorm(params_est["R0_est"], -0.25, 0.7, log = TRUE), # log(R0)
            dnorm(-params_est["dur_latent_est"], 0, 0.22, log = TRUE), # -log(dur_latent)
            dnorm(-params_est["dur_early_est"], 0, 0.22, log = TRUE), # -log(dur_early)
            dnorm(-params_est["dur_progress_est"], 0, 0.22, log = TRUE), # -log(dur_progress)
            dbeta(expit(params_est["ifr_est"]), 1.5, 200, log = TRUE) + params_est["ifr_est"] - 2 * log(exp(params_est["ifr_est"]) + 1), # logit(ifr)
            dbeta(expit(params_est["rho_death_est"]), 8, 2, log = TRUE) + params_est["rho_death_est"] - 2 * log(exp(params_est["rho_death_est"]) + 1) , # logit(rho_death)
            dexp(exp(params_est["phi_death_est"]), 1, log = TRUE) +  params_est["phi_death_est"], # -0.5 * log(phi_death)
            dnorm(exp(params_est["alpha0_est"]), 4, 2, log = TRUE) - log(2) + params_est["alpha0_est"], # log(alpha0)
            dbeta(expit(params_est["alpha1_est"]), 3, 1, log = TRUE) + params_est["alpha1_est"] - 2 * log(exp(params_est["alpha1_est"]) + 1), # logit(alpha1)
            dexp(exp(params_est["kappa_est"]), 1, log = T) +  params_est["kappa_est"]) # -0.5 * log(kappa)
      }

# return all three functions in a list
priors <- list(logprior = logprior,
               to_estimation_scale = to_estimation_scale,
               from_estimation_scale = from_estimation_scale)


#'
#' We now specify the MCMC transition kernel. In this simple example, we'll update
#' the model hyperparameters using a multivariate Metropolis algorithm. We'll tune
#' the algorithm using a global adaptive scheme (algorithm 4 in Andrieu and Thoms). We'll also initialize the parameters at random values, which is done by replacing the vector of parameters in the stem object with a function that returns a named vector of parameters.
#'
## ----mcmc_kern, echo = TRUE----------------------------------------------
# specify the initial proposal covariance matrix with row and column names
# corresponding to parameters on their estimation scales

par_initializer = function() {
  priors$from_estimation_scale(priors$to_estimation_scale(parameters) + rnorm(5, 0, 0.1))
}

# specify the kernel
mcmc_kern <-
        mcmc_kernel(
          parameter_blocks =
              list(parblock(
                  pars_nat = c("beta", "gamma", "nu_early", "mu_rec", "mu_death", "rho_death", "phi_death", "alpha0", "alpha1", "kappa"),
                  # par_est should have different names from pars_nat
                  pars_est = c("R0_est", "dur_latent_est", "dur_early_est", "dur_progress_est", "ifr_est", "rho_death_est", "phi_death_est", "alpha0_est", "alpha1_est", "kappa_est"),
                  priors = priors,
                  # alg = "mvnss",
                  alg = "mvnmh",
                  sigma = diag(0.01, 10),
                  initializer = par_initializer,
                  control =
                    # mvnss_control(stop_adaptation = 1e2))),
                    mvnmh_control(stop_adaptation = 1e4,
                                  scale_cooling = 0.8,
                                  scale_constant = 1,
                                  step_size = 0.4))),
          lna_ess_control = lna_control(bracket_update_iter = 5e3,
                                        joint_initdist_update = FALSE))

#'
#' We now run the MCMC algorithm to fit the model via ODEs. To perform inference with the LNA, simply change the `method` argument to `method="lna"`.
#'
## ----fit_mod, echo = TRUE------------------------------------------------

res <-
    fit_stem(stem_object = stem_object,
             method = "ode",
             mcmc_kern = mcmc_kern,
             iterations = 5e5,
             thinning_interval = 1000,
             print_progress = 1e4)

write_rds(res, "~/Desktop/res.rds")

#'
#' The `fit_stem` function returns a list with posterior samples, latent epidemic paths, and MCMC tuning parameters (e.g., global scaling parameter adapted in the MCMC). These can be accessed as follows:
#'
## ----access_res, echo = TRUE---------------------------------------------
runtime = res$results$runtime
mcmc_samples = res$results$posterior
ode_paths = res$results$ode_paths # or res$results$lna_paths if using the LNA for inference

#' # References
#'
#' Andrieu, C., and Thoms, J.. "A tutorial on adaptive MCMC." __Statistics and Computing__ 18.4 (2008): 343-373.
#'
#' Fintzi, J., Wakefield, J., & Minin, V. N. (2020). __A linear noise approximation for stochastic epidemic models fit to partially observed incidence counts.__ arXiv preprint arXiv:2001.05099.
#'
#' Gillespie, D. T. (1976). __A general method for numerically simulating the stochastic time__
#' __evolution of coupled chemical reactions.__ Journal of Computational Physics 22, 403–434.
#'
#' Murray, I., Adams, R. P., and MacKay, D. J. C. (2010). __Elliptical slice sampling.__ JMLR:
#' W&CP 9, 541–548.
