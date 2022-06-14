#A function that creates mutation events, necessary for change in Black Queen role of each species over 2000 generations
#Rates of mutation events: 
#Note! The mutation rate is controlled by the if(myevent < x ...) statement!  

#Realistic rate, from Cooper and Lenski, 2000: loss of 4.5e-03 catabolic functions per generation, measured over 2000 generations
#High rate, an order of magnitude above Cooper and Lenski: 0.045
#Low rate, an order of magnitude below Cooper and Lenski: 4.5e-04
#No mutations: 0

mutationevents <- function(mutationtimes, MutationRate){
  for(y in 2:2001){
    events <- runif(n = 20, min = 0, max = 1)
    for(i in 1:20){
      myevent <- events[i]
      if(myevent < MutationRate & mutationtimes[i] == 0){
        mutationtimes[i] <- round(y, -1)
      }
    }
  }
  return(mutationtimes)
}