run_linear_plots_meta_analysis_interaction <- function(exposure = "x", outcome = "y", seed = 124, replicates = 10, u = FALSE, n = 25000) {




run_internal_replicates <- function(replicates) {
set.seed = rnorm(replicates) + seed
value = rnorm(1, mean = 0, sd = 0.3)
linear = rbinom(1,1, 0.5)

mydata <- sim_mydata_outcomes_ma_interaction(n = n, 
                            linear = linear,
                            value = value)


z <- generate_all_sumstats_ma(mydata, exposure = exposure, outcome = outcome) %>%
  mutate(true = value)
z %>%
  mutate(exposure = exposure, outcome = outcome, linear = linear)
}


map_dfr(1:replicates, run_internal_replicates)



}




