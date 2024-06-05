sim_mydata_outcomes_ma_vQTL <- function(n = 25000,
                                          seed = 123,
                                          linear,
                                          
                                          value = value
                                          
){ 
  
  ########################
  ## Simulate the data 
  ########################
  set.seed(seed)
  ########################
  ## Simulate the data 
  ########################
  
  
  
  
    

    
    ########################
    ## Simulate the data 
    ########################
  
  
  if (linear == 0) {
    mydata <- tibble(
      g = rnorm(n),
      u = rnorm(n),
      g_rank = rank(g) / max(rank(g)),
      x = rnorm(n, mean = 0, sd = 1 + 0.3*g_rank) +  0.3 * g + 0.5 * u,
      y = rnorm(n) + value * 0.3*x^2 + u * 0.4) %>%
      as.data.frame()
    
  } else {
    mydata <- tibble(
      g = rnorm(n),
      u = rnorm(n),
      g_rank = rank(g) / max(rank(g)),
      x = rnorm(n, mean = 0, sd = 1 + 0.3*g_rank) +  0.3 * g + 0.5 * u,
      y = rnorm(n) + value * x + u * 0.4) %>%
      as.data.frame()
    
  }
  
  
  ########################
  ## return the sims
  ########################
  return(mydata)
}





