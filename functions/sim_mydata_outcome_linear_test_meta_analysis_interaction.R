sim_mydata_outcomes_ma_interaction <- function(n = 25000,
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
      x = rnorm(n) +  0.1 * g + 0.1 * u + -0.1*g*u,
      y = rnorm(n) + value * 0.3*x^2 + u * 0.4) %>%
      as.data.frame()
    
  } else {
    mydata <- tibble(
      g = rnorm(n),
      u = rnorm(n),
      x = rnorm(n) +  0.1* g + 0.1 * u + -0.1 * g * u,
      y = rnorm(n) + value * x + u * 0.4) %>%
      as.data.frame()
    
  }
  
  
  ########################
  ## return the sims
  ########################
  return(mydata)
}





