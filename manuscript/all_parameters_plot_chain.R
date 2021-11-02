tmp <-
  c("final_models/results/stemr_all_params/simulated_data_stemr_all_params.rds" = path("final_models/results/stemr_all_params/simulated_data_stemr_all_params.rds")) %>%
  map(~read_rds(.) %>%
        filter(time == 0) %>%
        select(starts_with("."), gamma, dur_infec, rho_death, phi_death, kappa,
               sigma_R0, sigma_ifr, sigma_alpha,
               S, E, I, R, D,
               R0_init, alpha_init, ifr_init)) %>%
  `names<-`(., names(.) %>% path_file() %>% path_ext_remove() %>% str_sub(end = -18)) %>%
  `names<-`(., model_name_conversion[names(.)]) %>%
  imap_dfr(~mutate(.x, model = .y)) %>%
  transmute(.chain = .chain,
            model = model,
            "$S_0$" = S / (S + E + I + R + D),
            "$\\tilde{I}_0$" = I / (E + I + R + D),
            "$\\tilde{R}_0$" = R / (E + R + D),
            "$\\tilde{D}_0$" = D / (E + D),
            "$1 / \\gamma$" = 1 / gamma,
            "$1 / \\nu" = dur_infec,
            "$\\rho$" = rho_death,
            "$1 / \\sqrt{\\phi}$" = exp(-log(phi_death)/2),
            "$1 / \\sqrt{\\kappa}$" = exp(-log(kappa)/2),
            "$\\exp{\\tilde{R}_{01}}$" = R0_init,
            "$\\expit{\\tilde{\\eta}_1}$" = ifr_init,
            "$\\exp{\\tilde{\\alpha}_1}$" = alpha_init,
            "$\\sigma_{R_0}$" = sigma_R0,
            "$\\sigma_{IFR}$" = sigma_ifr,
            "$\\sigma_{\\alpha}$" = sigma_alpha) %>%
  pivot_longer(-c(model, .chain)) %>%
  mutate(Distribution = "Posterior") %>%
  mutate(name = as_factor(name)) %>%
  mutate(name = fct_relabel(name, ~as.character(TeX(.)))) %>%
  filter(model == "Simulated Data")


tmp %>%
  ggplot(aes(value)) +
  facet_grid(.chain ~ name, labeller = label_parsed, scales = "free") +
  stat_halfeye(normalize = "xy")

