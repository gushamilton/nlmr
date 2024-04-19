generate_outcome_plots_linear_replicates_confounder <- function(exposure, outcome, seed = "123", reps = 10,
                                                                bux = 0.3,
                                                                buy = 0.3,
                                                                bgux = -0.1) {
  

  replicate_internally <- function(n_rep) {
  seed = seed + n_rep
  dat <-  sim_mydata_outcomes_confounder(n = 100000, seed = seed,
                                         bux = bux,
                                         buy = buy,
                                         bgux = bgux)
  dat_with_both <- generate_all_sumstats(data = dat, exposure = exposure, outcome = outcome, k = 10)
  
  ## redefine data set
  
  
  ## define sumstats tibble
  sum_stats_dat = dat_with_both$ss
  
  sum_stats_dat %>%
    mutate(replicate = n_rep)
  }
  
  d <- map_dfr(1:reps, replicate_internally)
 
  dat <-  sim_mydata_outcomes_confounder(n = 100000, seed = seed)
  
  ## make figures


  make_figures_replicates(d,exposure,outcome,dat,reps, "figures/linear/", paste0("_linear_effect_confounder_bux",bux,"_buy",buy,"_bgux",bgux,"_")) 
}





