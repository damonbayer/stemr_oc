time_interval_in_days <- 7

age_grouped_dat <- neg_line_list %>%
  filter(test_result != "other") %>%
  filter(date >= "2020-03-30",
         date <= "2020-10-26") %>%
  drop_na() %>%
  arrange(age) %>%
  # mutate(age_group = case_when(
  #                              age < 25 ~ "< 25",
  #                              25 <= age & age <= 65 ~ "25 - 65",
  #                              age > 65 ~ "> 65")) %>%
  mutate(age_group = case_when(
    age < 25 ~ "< 25",
    25 <= age & age <= 75 ~ "25 - 75",
    age > 75 ~ "> 75")) %>%
  mutate(age_group = fct_inorder(age_group)) %>%
  select(-id)

grouped_plot <- age_grouped_dat %>%
  count(date, test_result, age_group) %>%
  group_by(date) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(date, prop, color = age_group, linetype = test_result)) +
  geom_smooth(se = F) +
  theme_cowplot()

grouped_plot + facet_grid(age_group ~ test_result, scales = "free_y")
grouped_plot + facet_grid(test_result ~ age_group, scales = "free_y")
grouped_plot + facet_wrap(. ~ test_result, scales = "free_y")
grouped_plot + facet_wrap(. ~ age_group, scales = "free_y")

lumped_dat <- age_grouped_dat %>%
  group_by(lump = as.integer(floor((max(date) - date) / time_interval_in_days))) %>%
  mutate(lump = max(date)) %>%
  count(lump, test_result, age_group) %>%
  mutate(prop = n / sum(n))

lumped_plot_n <- lumped_dat %>%
  ggplot(aes(lump, n, color = age_group, linetype = test_result)) +
  geom_line() +
  theme_cowplot() +
  ylab("Count") +
  scale_y_continuous(labels = scales::comma)

lumped_plot_prop <- lumped_dat %>%
  ggplot(aes(lump, prop, color = age_group, linetype = test_result)) +
  geom_line() +
  theme_cowplot() +
  ylab("Proportion of all tests") +
  scale_y_continuous(labels = scales::percent)

# lumped_plot_n + facet_grid(age_group ~ test_result, scales = "free_y")
# lumped_plot_n + facet_grid(test_result ~ age_group, scales = "free_y")
lumped_plot_n + facet_wrap(. ~ test_result, scales = "free_y")
save_plot(filename = "~/Desktop/lumped_plot_n.png",
          plot =  lumped_plot_n + facet_wrap(. ~ test_result, scales = "free_y") + theme(legend.position = "bottom"),
          ncol = 2, nrow = 1)
# lumped_plot_n + facet_wrap(. ~ age_group, scales = "free_y")


# lumped_plot_prop + facet_grid(age_group ~ test_result, scales = "free_y")
# lumped_plot_prop + facet_grid(test_result ~ age_group, scales = "free_y")
lumped_plot_prop + facet_wrap(. ~ test_result, scales = "free_y")
save_plot(filename = "~/Desktop/lumped_plot_prop.png",
          plot =  lumped_plot_prop + facet_wrap(. ~ test_result, scales = "free_y") + theme(legend.position = "bottom"),
          ncol = 2, nrow = 1)
# lumped_plot_prop + facet_wrap(. ~ age_group, scales = "free_y")

positivity_plot <- lumped_dat %>%
  select(-prop) %>%
  arrange(lump, age_group, test_result) %>%
  count(lump, age_group, test_result, wt = n) %>%
  group_by(lump, age_group) %>%
  add_count(wt = n) %>%
  mutate(prop = n / nn) %>%
  filter(test_result == "positive") %>%
  ggplot(aes(lump, prop, color = age_group, group = age_group)) +
  geom_line() +
  scale_y_continuous(labels = scales::percent, name = "Positivity")

positivity_plot

save_plot(filename = "~/Desktop/lumped_plot_pos.png", plot = positivity_plot)
