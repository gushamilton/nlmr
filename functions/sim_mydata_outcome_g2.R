sim_mydata_outcomes_g2 <- function(n = 100000,
                                           mu = c(0, 0, 0, 0, 0), 
                                           sd = c(1, 1, 1, 1, 1), 
                                           r = 0, 
                                           varnames = c("g", "x", "y", "u", "v"),
                                           bv = 0,
                                           g2 = 0,
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
    
    mutate(x = (range01(g*10))^3 *g2 + bv * v + rnorm(n)) %>%
    mutate(y = 0 * x + bv * v + rnorm(n))
  
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
  return(mydata)
}

s <- sim_mydata_outcomes_g2(seed = 203, n = 100000,g2 = 5)
s %>%
  mutate(x = exp(g + 0.2*rnorm(1e5))) %>%
  ggplot(aes(x = g, y = x)) +
  geom_point()

s %>%
  mutate(x = exp(g + 0.2*rnorm(1e5))) %>%
  mutate(strat = ntile(g, 10)) %>%
  group_by(strat) %>%
  summarise(d = sd(x))
