library(tidyverse)
library(fs)
library(lubridate)

result_file_key <-
  tibble(full_path = dir_ls(path("/dfs6", "pub", "bayerd"))) %>%
  filter(path_ext(full_path) == "rds") %>%
  mutate(file_name = full_path %>% path_file() %>% path_ext_remove()) %>%
  mutate(model_name = file_name %>%
           str_sub(start = 28, end = -23) %>%
           str_remove("_$"),
         finish_time = file_name %>%
           str_sub(start = -19) %>%
           ymd_hms(),
         seed = file_name %>%
           str_sub(start = -22, end = -21) %>%
           str_remove("_") %>%
           as.integer()) %>%
  select(-file_name)

model_names <- unique(result_file_key$model_name)

save_dir <- path(path("/dfs6", "pub", "bayerd", "combined_results"))
dir_create(save_dir)
# current_model_name <- model_names[1]
for (current_model_name in model_names[-(1:13)]) {
# for (current_model_name in model_names) {
  print(current_model_name)
  model_paths <-
    result_file_key %>%
    filter(model_name == current_model_name) %>%
    arrange(seed) %>%
    pull(full_path)

  multi_chain_stem_fit <- read_rds(model_paths[1])

  multi_chain_stem_fit$stem_fit_list <-
    model_paths %>%
    map(read_rds) %>%
    map("stem_fit_list") %>%
    unname()

  last_iteration <- length(multi_chain_stem_fit$stem_fit_list[[1]]$results$posterior$data_log_lik)

  for (. in seq_along(multi_chain_stem_fit$stem_fit_list)) {
    multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$data_log_lik <- multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$data_log_lik[-last_iteration]
    multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$params_log_prior<- multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$params_log_prior[-last_iteration]
    multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$parameter_samples_nat <- multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$parameter_samples_nat[-last_iteration,]
    multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$parameter_samples_est <- multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$parameter_samples_est[-last_iteration,]
    multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$latent_paths <- multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$latent_paths[,,-last_iteration]
    multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$initdist_log_lik <- multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$initdist_log_lik[-last_iteration]
    multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$initdist_samples <- multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$initdist_samples[-last_iteration,]
    multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$initdist_draws[[1]] <- multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$initdist_draws[[1]][,-last_iteration]
    multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$tparam_log_lik <- multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$tparam_log_lik[-last_iteration]
    multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$tparam_samples <- multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$tparam_samples[,,-last_iteration]
    multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$tparam_draws[[1]] <- multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$tparam_draws[[1]][,-last_iteration]
    multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$tparam_draws[[2]] <- multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$tparam_draws[[2]][,-last_iteration]
    multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$tparam_draws[[3]] <- multi_chain_stem_fit$stem_fit_list[[.]]$results$posterior$tparam_draws[[3]][,-last_iteration]
  }

  names(multi_chain_stem_fit$stem_fit_list) <- NULL

  multi_chain_stem_fit$n_chains <- length(multi_chain_stem_fit$stem_fit_list)
  multi_chain_stem_fit$n_iterations <- length(multi_chain_stem_fit$stem_fit_list[[1]]$results$posterior$data_log_lik)

  write_rds(multi_chain_stem_fit, path(save_dir, current_model_name, ext = "rds"))
}

