library(tidyverse)
library(fs)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source("~/Documents/uci_covid_modeling/code/SEIeIpRD/helper_functions.R")
source("~/Documents/uci_covid_modeling/code/SEIeIpRD/plot_functions.R")

model_objects <- read_rds("~/Documents/uci_covid_modeling/code/SEIeIpRD/oc/2020-08-12_2020-09-16/model_objects.rds")

init_states <- c(S = 3012564L, E = 14866L, Ie = 20907L, Ip = 20590L, R = 0L, D = 0L)
model_objects$popsize <- sum(init_states)

model_objects$frac_carrs <- 1 - init_states[["S"]] / model_objects$popsize
model_objects$frac_carrs_infec <- (init_states[["Ie"]] + init_states[["Ip"]]) / (init_states[["E"]] + init_states[["Ie"]] + init_states[["Ip"]])
model_objects$frac_infec_early <- init_states[["Ie"]] / (init_states[["Ie"]] + init_states[["Ip"]])

all.equal((1 - model_objects$frac_carrs) * model_objects$popsize, init_states[["S"]])
all.equal(model_objects$frac_carrs_infec * model_objects$frac_carrs * model_objects$popsize, (init_states[["Ie"]] + init_states[["Ip"]]))
all.equal(model_objects$frac_carrs_infec * model_objects$frac_carrs * model_objects$frac_infec_early * model_objects$popsize, init_states[["Ie"]])

write_rds(model_objects, "~/Documents/stemr_oc/fixed_init_sim/model_objects.rds")

control_list <- list(adapt_delta = 0.99,
                     max_treedepth = 20)
# Sample Posterior --------------------------------------------------------
oc_post <- stan(file = "~/Documents/stemr_oc/fixed_init_sim/SEIeIpRD.stan",
                data = model_objects,
                seed = 0,
                chains = 4,
                control = control_list)

write_rds(oc_post, "fixed_init_sim/oc_post.rds")

# Sample Prior ------------------------------------------------------------
model_objects_priors_only <- model_objects
model_objects_priors_only$priors_only <- T
model_objects_priors_only$forecast_in_days <- 0

oc_prior <- stan(file = "~/Documents/stemr_oc/fixed_init_sim/SEIeIpRD.stan",
                 data = model_objects_priors_only,
                 seed = 0,
                 chains = 4,
                 control = control_list)

write_rds(oc_prior, "fixed_init_sim/oc_prior_tmp.rds")
