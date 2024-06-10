generate_outcome_plots_g2<- function(exposure, outcome, seed = "123", reps = 10,
                                       g2 = 0.4,
                                     bux = 0,

                                                                buy = 0) {
  

  replicate_internally <- function(n_rep) {
  seed = seed + n_rep
  dat <-  sim_mydata_outcomes_g2(n = 100000, seed = seed,
                                   g2 = g2,
                                 buy = buy,

                                         bux = bux)
  dat_with_both <- generate_all_sumstats(data = dat, exposure = exposure, outcome = outcome, k = 10)
  
  ## redefine data set
  
  
  ## define sumstats tibble
  sum_stats_dat = dat_with_both$ss
  
  sum_stats_dat %>%
    mutate(replicate = n_rep)
  }
  
  d <- map_dfr(1:reps, replicate_internally)
 
  dat <-  sim_mydata_outcomes_g2(n = 100000, seed = seed, g2 = g2, bux = bux, buy = buy)
  
  ## make figures


  make_figures_replicates(d,exposure,outcome,dat,reps, "figures/linear/", paste0("exponential_bx=", g2, "_bux=", bux, "_buy=", buy)) 
}







