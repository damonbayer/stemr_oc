library(tidyverse)
library(lubridate)
library(zoo)
multi_chain_stem_fit <- read_rds("~/Documents/stemr_oc/final_models/results/2021-02-28.rds")

last_day <- ymd("2021-02-28")
my_seed <- 1
source("final_models/forecast_results/compile_model_for_forecast.R")

stem_object$dynamics$initializer[[1]]$fixed <- T
stem_object$dynamics$fixed_inits <- T

max_chain <- multi_chain_stem_fit$stem_fit_list %>%
  map("results") %>%
  map("posterior") %>%
  map("data_log_lik") %>%
  map_dbl(max) %>%
  which.max()

max_iteration <- which.max(multi_chain_stem_fit$stem_fit_list[[max_chain]]$results$posterior$data_log_lik)


simulation_parameters <- multi_chain_stem_fit$stem_fit_list[[max_chain]]$results$posterior$parameter_samples_nat[max_iteration,] %>%
  signif(2)

init_dist <- multi_chain_stem_fit$stem_fit_list[[max_chain]]$results$posterior$initdist_samples[max_iteration,]

tparam_draws <- list(
  multi_chain_stem_fit$stem_fit_list[[max_chain]]$results$posterior$tparam_draws[[1]][,max_iteration],
  multi_chain_stem_fit$stem_fit_list[[max_chain]]$results$posterior$tparam_draws[[2]][,max_iteration],
  multi_chain_stem_fit$stem_fit_list[[max_chain]]$results$posterior$tparam_draws[[3]][,max_iteration]) %>%
  map(~rollmean(., k = 5, fill = "extend"))

simulation_parameters <- c(R0_init = 1.2, gamma = 3.1, dur_infec = 0.68, rho_death = 0.92,
                           phi_death = 2800, alpha_init = 4.1, kappa = 1.9e+08, ifr_init = 0.0022,
                           sigma_R0 = 0.054, sigma_ifr = 0.086, sigma_alpha = 0.072)

init_dist <- c(S_0 = 3170103, I_0 = 3594, R_0 = 1, D_0 = 1, E_0 = 1993)


tparam_draws <- list(c(-0.0480564342976278, -0.0480564342976278, -0.0480564342976278,
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
                       -0.0556229033996398, -0.0556229033996398), c(-0.0858908049474203,
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
                       ), c(-0.914930643328379, -0.914930643328379, -0.914930643328379,
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
                            -0.71699072262182))

stem_object$dynamics$initdist_params <- stem_object$dynamics$initdist_priors <- stem_object$dynamics$initializer[[1]]$init_states <- stem_object$dynamics$initializer[[1]]$prior <- init_dist
set.seed(1)

simulation_parameters[["phi_death"]] <- exp(5.4)
simulation_parameters[["kappa"]] <- exp(7.0)

simulate_stem(stem_object = stem_object,
              nsim = 1,
              simulation_parameters = list(simulation_parameters),
              tparam_draws = list(tparam_draws),
              paths = F,
              method = "ode")$natural_paths[[1]][,c("S", "E", "I", "R", "D")]

simulate_stem(stem_object = stem_object,
              nsim = 1,
              simulation_parameters = list(simulation_parameters),
              tparam_draws = list(tparam_draws),
              paths = F,
              method = "ode")$paths[[1]][,c("E2I", "I2D")]

simulate_data <- function() {
  simulate_stem(stem_object = stem_object,
                nsim = 1,
                simulation_parameters = list(simulation_parameters),
                tparam_draws = list(tparam_draws),
                paths = F,
                method = "ode")$datasets[[1]]
}

tmp <-
  rerun(.n = 1000, simulate_data()) %>%
  map_dfr(~as_tibble(.) %>% mutate(seroprev_cases = ifelse(time == 20, seroprev_cases, NA))) %>%
  mutate(data = "simulated") %>%
  bind_rows(multi_chain_stem_fit$data %>%
              ungroup() %>%
              select(colnames(simulate_data())) %>%
              mutate(data = "real"))
tmp %>%
  pivot_longer(-c(time, data)) %>%
  ggplot(aes(time, value, color = data)) +
  facet_wrap(. ~ name, scales = "free_y") +
  geom_point()

set.seed(1)

simulated_data_raw <-
  simulate_data() %>%
  as_tibble() %>%
  mutate(seroprev_cases = ifelse(time == 20, seroprev_cases, NA))

multi_chain_stem_fit$data %>%
  mutate(cases = simulated_data_raw$cases,
         deaths = simulated_data_raw$deaths,
         seroprev_cases = simulated_data_raw$seroprev_cases) %>%
  write_csv("final_models/data/simulated_data.csv")
