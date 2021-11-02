library(tidyverse)
library(tidybayes)
library(scales)
library(cowplot)
library(fs)
library(stemr)
library(coda)
source('stemr_functions.R')
theme_set(cowplot::theme_minimal_hgrid())

# Plot Options
ci_width <- c(0.5, 0.8, 0.95)
date_breaks <- "2 month"
date_breaks3 <- "3 months"
date_labels <- "%b '%y"

model_name_conversion <- c(
  `2021-02-28` = "Orange County",
  double_initial_counts = "Sensitivity 1",
  higher_ifr = "Sensitivity 4",
  huntington_beach = "Huntington Beach",
  irvine = "Irvine",
  longer_infectious_period = "Sensitivity 2",
  lower_R0 = "Sensitivity 3",
  lower_alpha = "Sensitivity 5",
  santa_ana = "Santa Ana",
  simulated_data = "Simulated Data"
)


simulation_parameters <-
  c(
    R0_init = 1.2, gamma = 3.1, dur_infec = 0.68, rho_death = 0.92,
    phi_death = 2800, alpha_init = 4.1, kappa = 1.9e+08, ifr_init = 0.0022,
    sigma_R0 = 0.054, sigma_ifr = 0.086, sigma_alpha = 0.072
  )

tparam_draws <- list(c(
  -0.0480564342976278, -0.0480564342976278, -0.0480564342976278,
  -0.398021744475401, -0.384896089425589, -0.543795343827761, -0.582078110523899,
  0.0301231045564148, 0.488398585798905, 0.396747079866406, 0.420884755511259,
  0.383891086193548, -0.22091872367557, -0.314654206817478, -0.176784511029037,
  -0.205655309632317, -0.53631163504416, -0.00941627148562795,
  -0.396257074667674, -0.366828508772458, 0.0814377349838548, 0.405639795991884,
  0.00719164917412897, 0.238473292426527, 0.279098225661917, 0.0688741073829157,
  0.155193653405115, 0.572995325154714, 1.05146658672034, 1.17444913094143,
  1.45156399805211, 1.47236746207627, 1.31709382109943, 0.897423964452343,
  0.608559172144112, 0.194605949173914, -0.17847480291358, -0.959218117681321,
  -0.861251493815745, -0.781277110440929, -0.662928169141229, -0.437721869342718,
  0.363515160229099, 0.320598486576097, 0.44192173566626, -0.0556229033996398,
  -0.0556229033996398, -0.0556229033996398
), c(
  -0.0858908049474203,
  -0.0858908049474203, -0.0858908049474203, 0.549281683468169,
  0.64095571650229, 0.781765229965117, 0.769176425434089, 0.932088093586281,
  0.348909548357974, 0.372392244676009, 0.239760948631366, -0.0960032346428879,
  0.046240873067556, 0.458594132477201, 0.39592704468361, 0.665179480177885,
  0.5514788914257, -0.0824944240374928, -0.316211208023479, -0.404920036517901,
  -1.21457948915544, -0.277941394177284, 0.039080213094182, -0.0307054631042438,
  -0.33347693136861, -0.201565423463202, -1.10237442451264, -1.25931266996931,
  -0.871863709338789, -0.522727383532642, -0.354170498775835, 0.106893497890586,
  0.322197109852489, 0.389148605664507, 0.519343490213794, 1.11515536682002,
  0.839203825375812, 0.977322148327585, 0.889349662380904, 0.805555394206076,
  0.329363616142376, 0.568657193821404, 0.476059437578614, 0.0271112484047705,
  -0.0444357383554626, 0.218095975262093, 0.218095975262093, 0.218095975262093
), c(
  -0.914930643328379, -0.914930643328379, -0.914930643328379,
  -0.875687956084769, -0.98154886941342, -1.0334890479356, -0.764044118344194,
  -0.563683220759433, -0.475509270757854, 0.301666290834711, 0.73022926272768,
  0.611410649413829, 0.47888514306843, 0.215872761764973, -0.650744044711817,
  -1.16139063773596, -1.29025693931278, -1.22201695201995, -1.07762185420384,
  -0.82367033628249, -0.369145636187704, -0.332614570267693, -0.217258403909708,
  0.148356743916844, 0.315598140331652, -0.0134733875854982, 0.371685526276042,
  0.583129606723513, 0.241908205497748, 0.476986622395759, 0.737287776508774,
  0.568368166379149, 0.322547648076811, 0.248470659746473, -0.239017616067521,
  -0.652609099352945, -0.981556477449624, -0.727415478310115, -0.760055061611517,
  -0.910689845806744, -0.585595352768983, -0.292685273036282, -0.490805882973428,
  -0.744993269412262, -0.707686497904077, -0.71699072262182, -0.71699072262182,
  -0.71699072262182
))


stemr_all_time_varying_params <-
  dir_ls("final_models/results/stemr_all_params/") %>%
  tail(2) %>%
  map(~read_rds(.) %>%
        select(starts_with("."), date, beta_t, dur_infec, S, E, I, D, ifr_t, alpha0_t, E2I, cases, deaths)) %>%
  `names<-`(., names(.) %>% path_file() %>% path_ext_remove() %>% str_sub(end = -18)) %>%
  `names<-`(., model_name_conversion[names(.)]) %>%
  imap_dfr(~mutate(.x, model = .y)) %>%
  mutate(popsize = case_when(model == "Irvine" ~ 264975L,
                             model == "Santa Ana" ~ 198471L,
                             model == "Huntington Beach" ~ 358874L,
                             TRUE ~ 3175692L),
         initial_cases = case_when(model == "Irvine" ~ 55L,
                                   model == "Santa Ana" ~ 56L,
                                   model == "Huntington Beach" ~ 35L,
                                   TRUE ~ 760L)) %>%
  mutate(log_popsize = log(popsize)) %>%
  replace_na(list(cases = 0,
                  deaths = 0)) %>%
  group_by(.chain, .iteration, .draw, model) %>%
  mutate(cumulative_cases = cumsum(cases) + initial_cases,
         cumulative_deaths = cumsum(deaths)) %>%
  ungroup() %>%
  transmute(
    date = date,
    model = model,
    R0 = exp(log(beta_t) + log_popsize + log(dur_infec)),
    Rt = R0 * S / popsize,
    ifr_t = ifr_t,
    alpha0_t = alpha0_t,
    latent_case_ratio = E2I / cases,
    cumulative_incidence = popsize - S,
    cumalitive_latent_case_ratio = cumulative_incidence / cumulative_cases,
    prevalence = E + I,
    cumulative_latent_deaths = D,
    cumulative_cases,
    cumulative_deaths
  ) %>%
  filter(model == "Simulated Data")


stemr_all_time_varying_params %>%
  select(date) %>%
  distinct() %>%
  arrange(date) %>%
  mutate(R0 = exp(log(simulation_parameters[["R0_init"]]) + cumsum(c(0, tparam_draws[[1]])) * simulation_parameters[["sigma_R0"]]))



# R0 Plots ----------------------------------------------------------------
R0_plot <-
  stemr_all_time_varying_params %>%
  select(R0, date) %>%
  group_by(date) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  ggplot(aes(date, R0, ymin = .lower, ymax = .upper)) +
  # facet_wrap(. ~ model) +
  # geom_interval(
  #   show_point = T,
  #   shape = "_",
  #   point_color = "black",
  #   fatten_point = 1
  # ) +
  geom_lineribbon(
    color = "steelblue4") +
  geom_line(data = stemr_all_time_varying_params %>%
              select(date) %>%
              distinct() %>%
              arrange(date) %>%
              mutate(R0 = exp(log(simulation_parameters[["R0_init"]]) + cumsum(c(0, tparam_draws[[1]])) * simulation_parameters[["sigma_R0"]])),
            mapping = aes(date, R0, ymax = NULL, ymin = NULL),
            linetype = "dotted",
            size = 1.5) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_brewer() +
  scale_color_brewer() +
  scale_y_continuous(
    name = expression(R[0]),
    breaks = seq(0, 3, by = 0.25),
    limits = c(.75, 2.25)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks3,
               date_labels = date_labels) +
  theme(legend.position = "none") +
  ggtitle("Posterior Basic Reproductive Number",
          subtitle = "Dotted line indicates true value used in simulation")



# Rt Plot -----------------------------------------------------------------
true_Rt <-
  structure(list(date = structure(c(
    18350, 18357, 18364, 18371,
    18378, 18385, 18392, 18399, 18406, 18413, 18420, 18427, 18434,
    18441, 18448, 18455, 18462, 18469, 18476, 18483, 18490, 18497,
    18504, 18511, 18518, 18525, 18532, 18539, 18546, 18553, 18560,
    18567, 18574, 18581, 18588, 18595, 18602, 18609, 18616, 18623,
    18630, 18637, 18644, 18651, 18658, 18665, 18672, 18679, 18686
  ), class = "Date"), Rt = c(
    1.1970934208985, 1.19034141809117,
    1.18290128473596, 1.17467774661377, 1.14375467997919, 1.11360695246981,
    1.0742737699564, 1.03372294982667, 1.02794322587505, 1.04756696845126,
    1.06183013708595, 1.07709111368298, 1.08958742724442, 1.06583443963243,
    1.03666172759264, 1.01544717395473, 0.992944922377094, 0.953842089405216,
    0.94327242482438, 0.914044063606773, 0.887920928510153, 0.884586688290128,
    0.897616115227113, 0.892038735126295, 0.89825534577784, 0.907003542189562,
    0.90591294249381, 0.90945338450842, 0.934194376325279, 0.984942397583886,
    1.04535168634329, 1.12590477410074, 1.21324459940231, 1.29488065817964,
    1.34826702988543, 1.37769105544001, 1.36996639108429, 1.32619675081128,
    1.2207932977576, 1.12186504358994, 1.03062678218039, 0.95170315972873,
    0.891437745558806, 0.875647874107416, 0.861527338094747, 0.856693502858208,
    0.832356521759992, 0.811986215212738, 0.79513233476528
  )), row.names = c(
    NA,
    -49L
  ), class = c("tbl_df", "tbl", "data.frame"))


Rt_plot <-
  stemr_all_time_varying_params %>%
  select(Rt, date) %>%
  group_by(date) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  ggplot(aes(date, Rt, ymin = .lower, ymax = .upper)) +
  # facet_wrap(. ~ model) +
  # geom_interval(
  #   show_point = T,
  #   shape = "_",
  #   point_color = "black",
  #   fatten_point = 1
  # ) +
  geom_lineribbon(
    color = "steelblue4") +
  geom_line(data = true_Rt,
            mapping = aes(date, Rt, ymax = NULL, ymin = NULL),
            linetype = "dotted",
            size = 1.5) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_brewer() +
  scale_color_brewer() +
  scale_y_continuous(
    name = expression(R[t]),
    breaks = seq(0, 3, by = 0.25),
    limits = c(.50, 1.75)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks3,
               date_labels = date_labels) +
  theme(legend.position = "none") +
  ggtitle("Posterior Effective Reproductive Number",
          subtitle = "Dotted line indicates true value used in simulation")




# IFR Plots ---------------------------------------------------------------
ifr_t_plot <-
  stemr_all_time_varying_params %>%
  select(date, ifr_t) %>%
  group_by(date) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  ggplot(aes(date, ifr_t, ymin = .lower, ymax = .upper)) +
  geom_lineribbon(
    color = "steelblue4") +
  geom_line(data = stemr_all_time_varying_params %>%
              select(date) %>%
              distinct() %>%
              arrange(date) %>%
              mutate(ifr_t = expit(logit(simulation_parameters[["ifr_init"]]) + cumsum(c(0, tparam_draws[[2]])) * simulation_parameters[["sigma_ifr"]])),
            mapping = aes(date, ifr_t, ymax = NULL, ymin = NULL),
            linetype = "dotted",
            size = 1.5) +
  scale_fill_brewer() +
  # scale_color_brewer() +
  theme(legend.position = "none") +
  scale_y_continuous(
    name = "IFR",
    labels = function(.)
      scales::percent(., accuracy = 0.1),
    # trans = scales::logit_trans(),
    breaks = seq(0, 0.05, by = 0.001),
    limits = c(0, NA)
  ) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks3,
               date_labels = date_labels) +
ggtitle("Posterior Infection Fatality Ratio",
        subtitle = "Dotted line indicates true value used in simulation")





# alpha_t_plots -----------------------------------------------------------
alpha_t_plot <-
  stemr_all_time_varying_params %>%
  select(date, alpha0_t) %>%
  group_by(date) %>%
  median_qi(.width = ci_width) %>%
  drop_na() %>%
  ggplot(aes(date, alpha0_t, ymin = .lower, ymax = .upper)) +
  geom_lineribbon(
    color = "steelblue4") +
  geom_line(data = stemr_all_time_varying_params %>%
              select(date) %>%
              distinct() %>%
              arrange(date) %>%
              mutate(alpha0_t = exp(log(simulation_parameters[["alpha_init"]]) + cumsum(c(0, tparam_draws[[3]])) * simulation_parameters[["sigma_alpha"]])),
            mapping = aes(date, alpha0_t, ymax = NULL, ymin = NULL),
            linetype = "dotted",
            size = 1.5) +
  scale_fill_brewer() +
  # scale_color_brewer() +
  theme(legend.position = "none") +
  scale_y_continuous(name = expression(alpha),
                     # scale_y_continuous(name = expression(alpha[0]),
                     # trans = scales::log_trans(),
                     breaks = 0:5) +
  scale_x_date(name = "Date",
               date_breaks = date_breaks3,
               date_labels = date_labels) +
  ggtitle("Posterior \u03B1",
        subtitle = "Dotted line indicates true value used in simulation")



# Save Plots --------------------------------------------------------------
main_posterior_results_plot <-
  plot_grid(
    R0_plot + theme(legend.position = c(0.05, 0.825)) + scale_fill_brewer(name = "Credibility", labels = function(.) percent(as.numeric(.))),
    Rt_plot,
    ifr_t_plot,
    alpha_t_plot,
    align = "hv",
    nrow = 2,
    ncol = 2
  )

save_plot(
  filename = "~/Documents/oc_covid19_stemr_manuscript/figures/simulated_data_plots_time_varying.pdf",
  plot = main_posterior_results_plot,
  ncol = 2,
  nrow = 2,
  device = cairo_pdf
)

