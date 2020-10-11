res <- read_rds("multi_chain_stem_fit.rds")


res[[1]]
all.equal(res[[1]]$dynamics, res[[2]]$dynamics)
all.equal(res[[1]]$dynamics$parameters, res[[2]]$dynamics$parameters)
all.equal(res[[1]]$dynamics$stoich_matrix_ode, res[[2]]$dynamics$stoich_matrix_ode)


all.equal(res[[1]]$measurement_process, res[[2]]$measurement_process)
all.equal(res[[1]]$measurement_process$censusmat, res[[2]]$measurement_process$censusmat)


all.equal(res[[1]]$results, res[[2]]$results)

all.equal(res[[1]]$restart, res[[2]]$restart)

length(unique(res))
