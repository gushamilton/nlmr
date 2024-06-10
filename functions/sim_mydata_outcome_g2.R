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
  
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  
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
    mutate(y = 0 * x + buy * v + rnorm(n)) 
  
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

sim_mydata_outcomes_g2(seed = 123, bx = 1, exp =3) %>%
  ggplot(aes(x = g, y = x)) +
  geom_point()






sim_mydata_outcomes_g2(seed = 1333, g2 = 0.1, bux = -0.4, buy = 0, n = 1e5) %>%
  mutate(strata = ntile(g, 10000)) %>%
  mutate(dr = as.numeric(generate_ranked_strata(g,x,10))) %>%
  filter(dr == 1) %>%
  ggplot(aes(x = g, y = x, colour = strata)) +
  geom_point() +
  geom_smooth(method = "lm")

s <- sim_mydata_outcomes_g2(seed = 123, g2 = 0.05, bux = 0.4) 

glance(lm(x ~g, data = s))
sim_mydata_outcomes_g2(seed = 123, g2 = 0.5, bux = 0.4) %>%
  mutate(strat = ntile(g, 10)) %>%
  group_by(strat) %>%
  summarise(g = sd(x))

