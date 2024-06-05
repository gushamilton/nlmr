generate_outcome_plots_g2<- function(exposure, outcome, seed = "123", reps = 10,
                                       g2 = 0.4,

                                                                bv = 1) {
  

  replicate_internally <- function(n_rep) {
  seed = seed + n_rep
  dat <-  sim_mydata_outcomes_g2(n = 100000, seed = seed,
                                   g2 = g2,

                                         bv = bv)
  dat_with_both <- generate_all_sumstats(data = dat, exposure = exposure, outcome = outcome, k = 10)
  
  ## redefine data set
  
  
  ## define sumstats tibble
  sum_stats_dat = dat_with_both$ss
  
  sum_stats_dat %>%
    mutate(replicate = n_rep)
  }
  
  d <- map_dfr(1:reps, replicate_internally)
 
  dat <-  sim_mydata_outcomes_g2(n = 100000, seed = seed, g2 = g2, bv = bv)
  
  ## make figures


  make_figures_replicates(d,exposure,outcome,dat,reps, "figures/linear/", paste0("_linear_effect_g2_", bv, "_bv", g2, "_g2")) 
}



s <- generate_outcome_plots_g2("x", "y", reps = 40, g2 = 5, bv = 1, seed = 123)

lm()



