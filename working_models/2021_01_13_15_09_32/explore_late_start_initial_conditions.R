library(lubridate)
library(scales)
popsize <- 3175692L

target_date <- ymd("2020-08-16")

target_E_I <- oc_data %>% filter(date >= target_date - weeks(4), date <= target_date) %>% mutate(cases = cases * 7) %>% pull(cases) %>% sum()
target_E <- round(target_E_I / 4)
target_I <- target_E_I - target_E
target_R <- round(0.115 * popsize)
target_D <- oc_data %>% filter(date <= target_date) %>% pull(deaths) %>% sum()
target_S <- popsize - (target_E + target_I + target_R + target_D)


init_infected <- c(S = target_S, E = target_E, I = target_I, R = target_R, D = target_D)
C <- 1
n <- 1e5

my_function <- function(C) norm(as.matrix(quantile(x = extraDistr::rdirmnom(n = n, size = popsize, alpha = init_infected / C)[,4], probs = c(0.025, 0.975)) - popsize * c(0.105, 0.124)), type = "F")

tibble(x = 10^(seq(1, 6, by = 0.5))) %>%
  mutate(y = map_dbl(x, my_function)) %>%
  arrange(y)

res <- optimize(f = my_function, interval = c(100, 1000))

sim_data <- extraDistr::rdirmnom(n = n, size = popsize, alpha = init_infected / res$minimum) %>%
  `colnames<-`(names(init_infected)) %>%
  as_tibble() %>%
  pivot_longer(everything()) %>%
  mutate(name = fct_inorder(name))

sim_data %>%
  ggplot(aes(value, group = name, fill = name)) +
  facet_wrap(. ~ name, scales = "free") +
  geom_density() +
  cowplot::theme_cowplot() +
  scale_x_continuous(labels = comma) +
  geom_vline(data = enframe(init_infected) %>%
               add_row(tibble(
                 name = c("R", "R"),
                 value =  popsize * c(0.105, 0.124))),
             mapping = aes(xintercept = value)) +
    geom_vline(data = sim_data %>%
                 group_by(name) %>%
                 summarize(quant = quantile(value,  probs = c(0.025, 0.975))),
               mapping = aes(xintercept = quant), color = "red")



