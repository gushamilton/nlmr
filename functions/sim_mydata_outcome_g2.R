sim_mydata_outcomes_g2 <- function(n = 100000,
                                           mu = c(2, 0, 0, 0, 0), 
                                           sd = c(1, 1, 1, 1, 1), 
                                           r = 0, 
                                           varnames = c("g", "x", "y", "u", "v"),
                                           buy = 0,
                                     bux = 0,
                                           bx = 0,
                                   exp =1,
                                           seed = seed
){ 
  
  ########################
  ## Simulate the data 
  ########################
  set.seed(seed)
  
  
  ########################
  ## Simulate the data 
  ########################
  mydata <- faux::rnorm_multi(n = n, 
                              mu = mu,
                              sd = sd,
                              r = r, 
                              varnames = varnames,
                              empirical = FALSE) %>%
    
    mutate(x = g^exp *bx + bux * v + rnorm(n)) %>%
    mutate(y = 0.2 * x + buy * v + rnorm(n)) 

  ## return the sims
  ########################
  return(mydata)
}


