sim_mydata_outcomes_vQTL <- function(n = 100000,
                                           mu = c(0, 0, 0, 0, 0), 
                                           sd = c(1, 1, 1, 1, 1), 
                                           r = 0, 
                                           varnames = c("g", "x", "y", "u", "v"),
                                           bux = 0,
                                           bx = 0,
                                            buy = 0,
                                           g_sd = 0,
                                     u_sd = 0,
                                           seed = seed
){ 
  
  ########################
  ## Simulate the data 
  ########################
  set.seed(seed)
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  ########################
  ## Simulate the data 
  ########################
  mydata <- faux::rnorm_multi(n = 100000, 
                              mu = mu,
                              sd = sd,
                              r = r, 
                              varnames = varnames,
                              empirical = FALSE) %>%
    mutate(g_rank = range01(g),
           u_rank = range01(u)) %>%
    mutate(x = bx * g + rnorm(n, mean = 0, sd =  u_rank*u_sd + g_rank *g_sd + 1) + bux * v) %>%
    mutate(y = 0 * x + buy * v + rnorm(n))
  
  ########################
  ## produce transformed 
  ## versions of the exposure data
 
  ########################
  ## produce transformed 
  ## versions of the outcome data
  ########################
  

  
  ########################
  ## return the sims
  ########################
  return(mydata)
}





