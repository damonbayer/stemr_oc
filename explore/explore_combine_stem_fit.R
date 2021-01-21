library(tidyverse)
multi_chain_stem_fit <- read_rds("fixed_init_sim/fixed_init_multi_chain_stem_fit_tmp.rds")


coerce_stem_fit_list <- function(stem_fit_list, param_blocks  = F ) {
  combined_stem_fit <- stem_fit_list[[1]]

  combined_stem_fit$dynamics$parameters <-
    stem_fit_list %>%
    map("dynamics") %>%
    map("parameters")

  combined_stem_fit$measurement_process$censusmat <-
    stem_fit_list %>%
    map("measurement_process") %>%
    map("censusmat")

  combined_stem_fit$results$runtime <-
    stem_fit_list %>%
    map("results") %>%
    map("runtime")


  combined_stem_fit$results$posterior$lp <-
    stem_fit_list %>%
    map("results") %>%
    map("posterior") %>%
    map(function(posterior) {
      cbind(data_log_lik = posterior[["data_log_lik"]],
            params_log_prior = posterior[["params_log_prior"]],
            params_log_posterior = posterior[["data_log_lik"]] + posterior[["params_log_prior"]])
    })

  combined_stem_fit$results$posterior$data_log_lik <- NULL
  combined_stem_fit$results$posterior$params_log_prior <- NULL

  combined_stem_fit$results$posterior$parameter_samples_nat <-
    stem_fit_list %>%
    map("results") %>%
    map("posterior") %>%
    map("parameter_samples_nat")

  combined_stem_fit$results$posterior$parameter_samples_est <-
    stem_fit_list %>%
    map("results") %>%
    map("posterior") %>%
    map("parameter_samples_est")

  combined_stem_fit$results$posterior$latent_paths <-
    stem_fit_list %>%
    map("results") %>%
    map("posterior") %>%
    map("latent_paths")

  combined_stem_fit$restart$path$latent_path <-
    stem_fit_list %>%
    map("restart") %>%
    map("path") %>%
    map("latent_path")

  combined_stem_fit$restart$path$data_log_lik <-
    stem_fit_list %>%
    map("restart") %>%
    map("path") %>%
    map("data_log_lik")

  if (param_blocks) {
    combined_stem_fit$restart$param_blocks <-
      stem_fit_list %>%
      map("restart") %>%
      map("param_blocks")
  }

  combined_stem_fit
}


tmp <- coerce_stem_fit_list(stem_fit_list = multi_chain_stem_fit$stem_fit_list, param_blocks = F)
# diffs <- all.equal(multi_chain_stem_fit$stem_fit_list[[1]], multi_chain_stem_fit$stem_fit_list[[2]])
#
# still need to explore
# all.equal(multi_chain_stem_fit$stem_fit_list[[1]]$restart$param_blocks, multi_chain_stem_fit$stem_fit_list[[2]]$restart$param_blocks)

