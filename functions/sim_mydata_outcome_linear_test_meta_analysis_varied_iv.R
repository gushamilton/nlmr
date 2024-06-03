sim_mydata_outcomes_ma_varied <- function(n = 25000,
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
  
  if (linear == 0) {
    mydata <- tibble(
      g = rnorm(n),
      x = rnorm(n) + if_else(g < 0, 0.2 * g, 0.5 *g^2 ),
      y = rnorm(n) + value* 0.3 *x^2) %>%
      as.data.frame()
    
  } else {
    mydata <- tibble(
      g = rnorm(n),
      x = rnorm(n) + if_else(g < 0, 0.2 * g, 0.5 *g^2 ),
      y = rnorm(n) + value * x) %>%
      as.data.frame()
    
  }
  

  ########################
  ## return the sims
  ########################
  return(mydata)
}



