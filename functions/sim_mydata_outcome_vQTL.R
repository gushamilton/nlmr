sim_mydata_outcomes_vQTL <- function(n = 100000,
                                           mu = c(0, 0, 0, 0, 0), 
                                           sd = c(1, 1, 1, 1, 1), 
                                           r = 0, 
                                           varnames = c("g", "x", "y", "u", "v"),
                                           bv = 0,
                                           bx = 0,
                                           g_sd = 0,
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
    mutate(g_rank = rank(g) / max(rank(g))) %>%
    mutate(x = bx * g + rnorm(n, mean = 0, sd =  g_rank *g_sd + 1) + bv * v) %>%
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

s <- sim_mydata_outcomes_vQTL(seed = 203, n = 100000)

s %>%
  as_tibble() %>%
  mutate(strat = ntile(g, 10)) %>%
  group_by(strat) %>%
  summarise(mean = mean(x), sd = sd(x))

model <- lm(x ~g, data = s)
tidy(model)
bptest(model)


s %>%
  as_tibble() %>%
  mutate(strat = ntile(g, 1000)) %>%
  mutate(dr = generate_ranked_strata(g,x,100)) %>%
  filter(dr == 100 | dr == 1) %>%
  # filter(strat == 1000 | strat == 1) %>%
  ggplot() +
  geom_smooth(aes(x = g, y = x, group = dr))  +
  geom_point(aes(x = g, y = x, colour = strat))

