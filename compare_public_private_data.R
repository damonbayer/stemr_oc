library(tidyverse)
library(lubridate)
library(cowplot)
public_data <- read_csv("/Users/damon/Desktop/public_OC_data.csv") %>%
  select(-county) %>%
  filter(date >="2020-03-30",
         date <= "2021-02-28") %>%
  group_by(time = as.numeric(floor((date - ymd("2020-03-30")) / 7 + 1))) %>%
  summarize(start_date = min(date),
            end_date = max(date),
            cases = sum(cases),
            tests = sum(tests),
            deaths = sum(deaths)) %>%
  mutate(source = "public")

private_data <- read_rds("~/Documents/stemr_oc/final_models/results/2021-02-28.rds")$data %>%
  mutate(source = "private") %>%
  select(names(public_data))



bind_rows(
  public_data,
  private_data) %>%
  mutate(positivity = cases / tests) %>%
  ggplot(aes(end_date, deaths, group = source, color = source)) +
  geom_point() +
  geom_line() +
  theme_minimal_grid()


bind_rows(
  public_data,
  private_data) %>%
  mutate(positivity = cases / tests) %>%
  ggplot(aes(end_date, cases, group = source, color = source)) +
  geom_point() +
  geom_line() +
  theme_minimal_grid()

bind_rows(
  public_data,
  private_data) %>%
  mutate(positivity = cases / tests) %>%
  ggplot(aes(end_date, positivity, group = source, color = source)) +
  geom_point() +
  geom_line() +
  theme_minimal_grid()
