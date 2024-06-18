generate_outcome_plots_g2<- function(exposure, outcome, seed = "123", reps = 10,
                                       bx = bx,
                                     exp = exp,
                                     bux = bux,
                                     n = n,

                                                                buy = buy) {
  

  replicate_internally <- function(n_rep) {
  seed = seed + n_rep
  dat <-  sim_mydata_outcomes_g2(n = n, seed = seed,
                                   bx = bx,
                                 exp = exp,
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
 
  dat <-  sim_mydata_outcomes_g2(n = n,seed = seed, bx = bx, bux = bux,exp = exp, buy = buy)
  
  ## make figures


  make_figures_replicates(d,exposure,outcome,dat,reps, "figures/linear/", paste0("exponential_bx=", bx, "_exp=",exp, "_bux=", bux, "_buy=", buy)) 
}


