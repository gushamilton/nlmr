sim_mydata_outcomes_confounder_RINT <- function(n = 100000,
                       mu = c(0, 0, 0, 0), 
                       sd = c(1, 1, 1, 1), 
                       r = c(0.3, 0.02, 0, 0.1, 0, 0), 
                       varnames = c("g", "x", "y", "u"),
                       bux = 0.3,
                       buy = 0.3,
                       bgux = -0.1,
                       seed = seed
                       ){ 
  
  ########################
  ## Simulate the data 
  ########################
  set.seed(seed)
  
  ########################
  ## Simulate the data 
  ########################
  mydata <- faux::rnorm_multi(n = 100000, 
                              mu = mu,
                              sd = sd,
                              r = r, 
                              varnames = varnames,
                              empirical = FALSE) %>%
    mutate(x = 0.3 * g + bux * u + bgux * g * u + rnorm(1e5)) %>%
    mutate(y = 0 * x + buy * u + rnorm(1e5))
  
  ########################
  ## produce transformed 
  ## versions of the exposure data
  ########################
  ## UNIFORM
  a = mydata$x
  a = (a - mean(a) )/ sd(a)
  mydata$x_uni = pnorm( a, mean = 0, sd = 1)
  
  ## GAMMA
  a = a + abs( min(a) )
  mydata$x_gamma = pgamma( a , shape = 3)
  
  ## EXPONENTIAL
  mydata$x_exp = pexp(mydata$x_uni, rate = 4)
  
  ## BETA
  a = a / max(a)
  mydata$x_beta = pbeta(a, shape1 = 4, shape2 = 3)
  
  ########################
  ## produce transformed 
  ## versions of the outcome data
  ########################
  
  
  ## UNIFORM
  a = mydata$y
  a = (a - mean(a) )/ sd(a)
  mydata$y_uni = pnorm( a, mean = 0, sd = 1)
  
  ## GAMMA
  a = a + abs( min(a) )
  mydata$y_gamma = pgamma( a , shape = 3)
  
  ## EXPONENTIAL
  mydata$y_exp = pexp(mydata$y_uni, rate = 4)
  
  ## BETA
  a = a / max(a)
  mydata$y_beta = pbeta(a, shape1 = 4, shape2 = 3)
  
  ########################
  ## return the sims
  ########################
  mydata <-mydata %>%
    mutate(across(contains("x"), rntransform))
  return(mydata)
}



