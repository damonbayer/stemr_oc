library(rstan)
files <- dir(system.file("include", "stan", "math", "prim",
                         package = "StanHeaders"),
             pattern = "hpp$", recursive = TRUE)
functions <- sub("\\.hpp$", "",
                 sort(unique(basename(files[dirname(files) != "."]))))

N <- model_objects$lumped_ochca_covid$new_tests[1]
x <- N / 2
alpha = 4.25
beta = 9.37

stanFunction("beta_binomial_lpmf", x = x, N = N, alpha = alpha, beta = beta) %>% sprintf("%.100f", .)
extraDistr::dbbinom(x = x, size = N, alpha = alpha, beta = beta, log = T) %>% sprintf("%.100f", .)
