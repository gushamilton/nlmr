sim_mydata_outcomes_final<- function(n = 100000,
                       seed = seed
                       ){ 
  
  ########################
  ## Simulate the data 
  ########################
  set.seed(seed)
  
  ########################
  ## Simulate the data 
  ########################
  mydata <- tibble(g = rnorm(n),
                   x = rnorm(n) + g * 0.3,
                   y = rnorm(n))
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
  

  mydata$y_uni = mydata$y + 0.3 * mydata$x_uni
  
  ## GAMMA
 
  mydata$y_gamma = mydata$y + 0.3 * mydata$x_gamma
  
  ## EXPONENTIAL
  mydata$y_exp = mydata$y + 0.3 * mydata$x_exp
  
  ## BETA
  mydata$y_beta = mydata$y + 0.3 * mydata$x_beta
  
  
  mydata$y = mydata$y + 0.3 * mydata$x
  
  mydata <- as.data.frame(mydata)
  ########################
  ## return the sims
  ########################
  return(mydata)
}

s <- sim_mydata_outcomes_final(100000, seed = 123)

s





bmi <- tibble( g = rnorm(10000),
weight = rnorm(10000, mean = 30, 3) ,
height = rnorm(10000, 180, 15) + g * 10,
bmi = weight/height^2)
ggplot(bmi) +
  geom_histogram(aes(x = bmi))


ggplot(bmi) +
  geom_point(aes(x =g, y = bmi)) +
  geom_smooth(aes(x =g, y = bmi))

                 