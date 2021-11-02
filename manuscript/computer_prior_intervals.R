# dnorm(params_est["R0_init_est"], 0, 0.225, log = TRUE), # log(R0)
# dnorm(-params_est["dur_latent_est"],  log(3/7), 0.25, log = TRUE), # -log(dur_latent)
# dnorm(params_est["dur_infec_est"], mean = log(7/7), sd = 0.225, log = T),
# dbeta(expit(params_est["rho_death_est"]), 200, 20, log = TRUE) + params_est["rho_death_est"] - 2 * log(exp(params_est["rho_death_est"]) + 1) , # logit(rho_death)
# dexp(exp(params_est["phi_death_est"]), 1, log = TRUE) +  params_est["phi_death_est"], # -0.5 * log(phi_death)
# dtnorm(exp(params_est["alpha0_init_est"]), mean = 4, sd = 0.5, a = 0, log = TRUE) + params_est["alpha0_init_est"], # log(alpha0)
# dexp(exp(params_est["kappa_est"]), 1, log = T) +  params_est["kappa_est"], # -0.5 * log(kappa)
# dbeta(expit(params_est["ifr_init_est"]), 30, 10000, log = TRUE) + params_est["ifr_init_est"] - 2 * log(exp(params_est["ifr_init_est"]) + 1), # logit(ifr_init)
# dtnorm(exp(params_est["sigma1_est"]), mean = 0.05, sd = 0.01, a = 0, log = TRUE) + params_est["sigma1_est"], # log(sigma1)
# dtnorm(exp(params_est["sigma2_est"]), mean = 0.08, sd = 0.01, a = 0, log = TRUE) + params_est["sigma2_est"], # log(sigma2)
# dtnorm(exp(params_est["sigma3_est"]), mean = 0.05, sd = 0.01, a = 0, log = TRUE) + params_est["sigma3_est"]  # log(sigma3)

library(extraDistr)

qs <- c(0.5, 0.05, 0.95)

library(glue)
library(tidyverse)
rq <- scales::comma(qtnorm(qs, mean = 2, sd = 0.5, a = 0))
glue("{rq[1]} ({rq[2]}, {rq[3]})")


rq <- scales::comma(qexp(qs, rate = 1))
glue("{rq[1]} ({rq[2]}, {rq[3]})")



rq <- scales::comma(qbeta(qs, 30, 10000))
glue("{rq[1]} ({rq[2]}, {rq[3]})")

rq <- scales::comma(qbeta(qs, 60, 9000))
glue("{rq[1]} ({rq[2]}, {rq[3]})")


rq <- scales::comma(exp(qnorm(qs, 0, 0.225)))
glue("{rq[1]} ({rq[2]}, {rq[3]})")

# scales::comma(c(log(3/7), 0.25), accuracy = 0.01)
scales::comma(c(log(0.5), 0.4), accuracy = 0.01)
rq <- scales::comma(exp(qnorm(qs, log(0.5), 0.4)))
glue("{rq[1]} ({rq[2]}, {rq[3]})")

scales::comma(c( 6.02067492520774, 0.549603697818656), accuracy = 0.01) %>% str_c(collapse = ", ") %>% cat()
rq <- scales::comma(stemr::expit(qnorm(qs, 6.02067492520774, 0.549603697818656)))
glue("{rq[1]} ({rq[2]}, {rq[3]})")

scales::comma(c(0.619153500651619, 0.0292277802326922), accuracy = 0.01) %>% str_c(collapse = ", ") %>% cat()
rq <- scales::comma(stemr::expit(qnorm(qs, 0.619153500651619, 0.0292277802326922)))
glue("{rq[1]} ({rq[2]}, {rq[3]})")

scales::scientific(c(-7.89766815072691, 9.8518767486955e-05)) %>% str_c(collapse = ", ") %>% cat()
rq <- scales::scientific(stemr::expit(qnorm(qs, -7.89766815072691, 9.8518767486955e-05)))
glue("{rq[1]} ({rq[2]}, {rq[3]})")

scales::scientific(c(-7.89729647259588, 9.84821737399636e-6), accuracy = 0.00001) %>% str_c(collapse = ", ") %>% cat()
rq <- scales::scientific(stemr::expit(qnorm(qs, -7.89729647259588, 9.84821737399636e-6)))
glue("{rq[1]} ({rq[2]}, {rq[3]})")


# -------------------------------------------------------------------------

scales::comma(c(5.3252274084003, 0.12208870486078), accuracy = 0.01) %>% str_c(collapse = ", ") %>% cat()
rq <- scales::comma(stemr::expit(qnorm(qs, 5.3252274084003, 0.12208870486078)))
glue("{rq[1]} ({rq[2]}, {rq[3]})")

scales::comma(c(0.619525040689332, 0.0641635068561201), accuracy = 0.01) %>% str_c(collapse = ", ") %>% cat()
rq <- scales::comma(stemr::expit(qnorm(qs, 0.619525040689332, 0.0641635068561201)))
glue("{rq[1]} ({rq[2]}, {rq[3]})")

scales::scientific(c(-8.59062950948942, 0.000237679033324254)) %>% str_c(collapse = ", ") %>% cat()
rq <- scales::scientific(stemr::expit(qnorm(qs, -8.59062950948942, 0.000237679033324254)))
glue("{rq[1]} ({rq[2]}, {rq[3]})")

scales::scientific(c(-8.59044365315583, 0.000237642407309249), accuracy = 0.00001) %>% str_c(collapse = ", ") %>% cat()
rq <- scales::scientific(stemr::expit(qnorm(qs, -8.59044365315583, 0.000237642407309249)))
glue("{rq[1]} ({rq[2]}, {rq[3]})")





