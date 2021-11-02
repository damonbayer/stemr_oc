library(tidyverse)
library(scales)
library(cowplot)
latent_incidence_denom <- rev(c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000))
latent_incidence_label <- str_c("1:", scales::comma(latent_incidence_denom))
latent_incidence <-  1 / latent_incidence_denom

latent_incidence_test_pos_plot <-
  expand_grid(alpha = seq(0, 9, 1), latent_incidence) %>%
  mutate(test_pos = stemr::expit(alpha + stemr::logit(latent_incidence))) %>%
  mutate(alpha = as_factor(alpha),
         latent_incidence = as_factor(latent_incidence)) %>%
  ggplot(aes(latent_incidence, test_pos, group = alpha, color = alpha)) +
  geom_line() +
  geom_point() +
  # scale_color_brewer(name = expression(alpha["\u2113"]), palette = "RdBu") +
  scale_color_brewer(name = expression(alpha), palette = "RdBu") +
  scale_x_discrete(name = "Latent Incidence (cases:population)", labels = latent_incidence_label) +
  scale_y_continuous(name = "Test Positivity Probability", labels = scales::percent) +
  cowplot::theme_minimal_hgrid()

save_plot(filename = "~/Documents/oc_covid19_stemr_manuscript/figures/latent_incidence_test_pos_plot.pdf",
          plot = latent_incidence_test_pos_plot, ncol = 1, nrow = 1, base_asp = 2.5, device = cairo_pdf)
