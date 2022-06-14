#Black Queen cellulose degradation simulations
#Authors: Damien Finn, damien.finn@thuenen.de
#         Mario App, mario.app@thuenen.de

rm(list=ls()) # clean up Environment

#Necessary packages
library(deSolve)
library(vegan)
#Set seed
set.seed(1234)
#set directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#To run the Black Queen simulations, set the following three parameters:

# 1) MutationRate. This defines the rate by which loss of function mutants appear in the simulations. 0 means no mutations. 0.045 is a relatively high rate.
#                  Realistic rate, from Cooper and Lenski, 2000: loss of 4.5e-03 catabolic functions per generation, measured over 2000 generations
#                  High rate, an order of magnitude above Cooper and Lenski: 0.045
#                  Low rate, an order of magnitude below Cooper and Lenski: 4.5e-04
#                  No mutations: 0
# 2) MonteCarlos.  This defines the number of Monte Carlo simulations are run.
# 3) BulkSoil.     If T, the Bulk Soil environment is simulated. If F, the Rhizosphere environment is simulated.


source("BlackQueenSim.R") # load the model
  
res <- BlackQueenSim(MutationRate = 0, MonteCarlos = 1, BulkSoil = T) # initialise the model

#The functional group outputs are as follows:
#CP = Cellulytic Prototrophs
#CA = Cellulytic Auxotrophs
#NCP = Noncellulytic Prototrophs
#NCA = Noncellulytic Auxotrophs