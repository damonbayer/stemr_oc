fixed_init_multi_chain_stem_fit <- readRDS("~/Documents/stemr_oc/fixed_init_sim/fixed_init_multi_chain_stem_fit.rds")
# oc_post <- readRDS("~/Documents/stemr_oc/fixed_init_sim/oc_post.rds") # Something wrong with this one
oc_post <- read_rds("/Users/damon/Documents/uci_covid_modeling/code/SEIeIpRD/oc/2020-08-12_2020-09-16/oc_post.rds")
model_objects <- readRDS("~/Documents/stemr_oc/fixed_init_sim/model_objects.rds")

source("~/Documents/uci_covid_modeling/code/SEIeIpRD/plot_functions.R")

oc_post_stemr_params <- map(fixed_init_multi_chain_stem_fit$stem_fit_list,
function(stem_fit) {
  stem_fit$results$posterior$parameter_samples_nat %>%
    split(., row(.)) %>%
    map(~`names<-`(.,c("beta", "gamma", "nu_early", "mu_rec", "mu_death", "rho_death",
                       "phi_death", "alpha0", "alpha1", "kappa")))
}) %>%
  unlist(recursive = F)


stemr_pp <-
  simulate_stem(stem_object = stem_object, method = "ode", nsim = 8000,
                simulation_parameters = oc_post_stemr_params) %>%
  `$`(dataset) %>%
  map_dfr(as_tibble)

prepped_death_stemr <- select(stemr_pp, deaths, time)


prepped_pos_stemr <- select(stemr_pp, cases, time) %>%
  left_join(model_objects$lumped_ochca_covid %>%
              mutate(time = as.numeric((end_date - min(end_date) + 3) / 7)) %>%
              select(time, tests = new_tests)) %>%
  mutate(percent_positive = cases / tests)


prepped_death_stan <- prep_death(model_objects, oc_post, 0.95) %>%
  map(~mutate(., time = as.numeric((t - min(t) + 3) / 7)) %>%
        select(-t))

prepped_pos_stan <- prep_pos(model_objects, oc_post, 0.95) %>%
  map(~mutate(., time = as.numeric((t - min(t) + 3) / 7)) %>%
        select(-t))





# death_plot <-
  ggplot() +
  geom_lineribbon(data = bind_rows(mutate(prepped_death_stan$posterior, model = "stan"),
                                   mutate(prepped_death_stemr, model = "stemr")) %>%
                    group_by(time, model) %>%
                    median_qi(.width = 0.95),
                  mapping = aes(time, deaths, ymin = .lower, ymax = .upper, fill = model),
                  size = 1.5, alpha = 0.5) +
  geom_point(data = prepped_death_stan$observed,
             mapping = aes(time, observed), size=2,
             fill = "black") +
  xlab("Date") +
  ylab("Number of Deaths")

  # oc_post <- readRDS("~/Documents/stemr_oc/fixed_init_sim/oc_post.rds") # Something wrong with this one
# TODO: positivity plot
