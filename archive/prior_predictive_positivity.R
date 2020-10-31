stemr_rds <- "priors_only_sim/priors_only_multi_chain_stem_fit.rds"
stan_rds <- "fixed_init_sim/oc_prior.rds"

stemr_prior_pred <- map(read_rds(stemr_rds)$stem_fit_list,
                        function(stem_fit) {
                          stem_fit$results$posterior$parameter_samples_nat %>%
                              split(., row(.)) %>%
                              map(~`names<-`(.,c("beta", "gamma", "nu_early", "mu_rec", "mu_death", "rho_death",
                                                 "phi_death", "alpha0", "alpha1", "kappa")))
                          }) %>%
    unlist(recursive = F) %>%
    simulate_stem(stem_object = stem_object, method = "ode", nsim = 8000,
                  simulation_parameters = .) %>%
    `$`(dataset) %>%
    map_dfr(as_tibble) %>%
    mutate(t = model_objects$dates_pp[time * 7 / 3]) %>%
    select(deaths, cases, t)


  prior_pred_pos <- prep_pos(model_objects = model_objects, oc_post = read_rds(stan_rds))

  prior_pred_pos$stemr <- stemr_prior_pred %>%
    select(t, cases) %>%
    left_join(select(model_objects$lumped_ochca_covid, t = end_date, tests = new_tests)) %>%
    mutate(stemr = cases / tests) %>%
    select(t, stemr) %>%
    group_by(t) %>%
    mutate(.draw = row_number())

  prior_pred_pos$stan <- prior_pred_pos$posterior %>% rename(stan = percent_positive) %>% filter(t %in% prior_pred_pos$stemr$t) %>% group_by(t) %>% mutate(.draw = row_number())
  prior_pred_pos$posterior <- NULL


  left_join(prior_pred_pos$stan, prior_pred_pos$stemr) %>%
    select(-.draw) %>%
    pivot_longer(-t) %>%
    ggplot(aes(logit(value), color = name, fill = name, group = name)) +
    stat_halfeye(alpha = 0.5) +
    facet_wrap(. ~ t) +
    ggtitle("Prior Predictive Positiviy Slices") +
    theme(legend.position = "bottom") +
    ylab(NULL) +
    xlab(NULL)
