BlackQueenSim <- function(MutationRate, MonteCarlos, BulkSoil){
  
  #Setup ODEs ----
  #A function to format the changes in cellulose, glucose, amino acids and species biomass in the Cellulose Degradation Simulation
  #This must be formatted as characters and then converted back to an expression with eval(parse(...)) when running as ODE
  source("SetUpODEs.R") 
  
  #Mutation event function----
  #A function that creates mutation events, necessary for change in Black Queen role of each species over 2000 generations
  #Rates of mutation events: 
  #Note! The mutation rate is controlled by the if(myevent < x ...) statement!  
  
  #Realistic rate, from Cooper and Lenski, 2000: loss of 4.5e-03 catabolic functions per generation, measured over 2000 generations
  #High rate, an order of magnitude above Cooper and Lenksi: 0.045
  #Low rate, an order of magnitude below Cooper and Lenski: 4.5e-04
  #No mutations: 0
  source("mutationevents.R")
  
  #A function to create new mutation profiles that build upon each other and keep track of mutation over time
  #Note that species with two public good functions are randomly mutated to only cellulase or amino acid production
  #Species that have only one public good function are mutated to Cheaters
  
  source("mutationsprofile.R")
  
  #Create a bunch of data frames to store various information 
  
  #Storing cumulative biomass (ng)
  CBdf <- 
    data.frame('Cumulative.Biomass' = rep(0, MonteCarlos)) 
  #Storing number of mutation events that arise
  mutationcounts <- 
    data.frame('Mutations' = rep(0, MonteCarlos)) 
  Scounts <- 
    data.frame('S1' = rep(0, MonteCarlos), 'S2' = rep(0, MonteCarlos),
               'S3' = rep(0, MonteCarlos), 'S4' = rep(0, MonteCarlos),
               'S5' = rep(0, MonteCarlos), 'S6' = rep(0, MonteCarlos),
               'S7' = rep(0, MonteCarlos), 'S8' = rep(0, MonteCarlos), 
               'S9' = rep(0, MonteCarlos), 'S10' = rep(0, MonteCarlos),
               'S11' = rep(0, MonteCarlos), 'S12' = rep(0, MonteCarlos), 
               'S13' = rep(0, MonteCarlos), 'S14' = rep(0, MonteCarlos), 
               'S15' = rep(0, MonteCarlos), 'S16' = rep(0, MonteCarlos), 
               'S17' = rep(0, MonteCarlos), 'S18' = rep(0, MonteCarlos), 
               'S19' = rep(0, MonteCarlos), 'S20' = rep(0, MonteCarlos)) 
  #Final distribution of species
  Cell.Prototrophs.B <- 
    data.frame('CP.Copios.Biom' = rep(0, MonteCarlos), 
               'CP.Oligos.Biom' = rep(0, MonteCarlos))
  NonCell.Prototrophs.B <- 
    data.frame('NCP.Copios.Biom' = rep(0, MonteCarlos), 
               'NCP.Oligos.Biom' = rep(0, MonteCarlos))
  Cell.Auxotrophs.B <- 
    data.frame('CA.Copios.Biom' = rep(0, MonteCarlos),
               'CA.Oligos.Biom' = rep(0, MonteCarlos))
  NonCell.Auxotrophs.B <- 
    data.frame('NCA.Copios.Biom' = rep(0, MonteCarlos),
               'NCA.Oligos.Biom' = rep(0, MonteCarlos)) #Biomass of the various species types
  Cell.Prototrophs.P <- 
    data.frame('CP.Prop' = rep(0, MonteCarlos),
               'CP.Copios.Prop' = rep(0, MonteCarlos),
               'CP.Oligos.Prop' = rep(0, MonteCarlos))
  NonCell.Prototrophs.P <- 
    data.frame('NCP.Prop' = rep(0, MonteCarlos),
               'NCP.Copios.Prop' = rep(0, MonteCarlos),
               'NCP.Oligos.Prop' = rep(0, MonteCarlos))
  Cell.Auxotrophs.P <- 
    data.frame('CA.Prop' = rep(0, MonteCarlos),
               'CA.Copios.Prop' = rep(0, MonteCarlos),
               'CA.Oligos.Prop' = rep(0, MonteCarlos))
  NonCell.Auxotrophs.P <- 
    data.frame('NCA.Prop' = rep(0, MonteCarlos),
               'NCA.Copios.Prop' = rep(0, MonteCarlos),
               'NCA.Oligos.Prop' = rep(0, MonteCarlos)) #Initial proportions of the various species types
  Shannondf <- 
    data.frame('Shannon' = rep(0, MonteCarlos)) # Storing alpha-diversity results of communities at final t
  Simpsondf <- 
    data.frame('Simpson' = rep(0, MonteCarlos)) # Also storing Simpson at final t
  
  #Running the Monte Carlo simulations
  for(y in 1:MonteCarlos){
    start.time <- Sys.time() #To keep track of computational time
    # Step 1: Set up dataframes etc.----
    #This controls whether the species produces Public Goods (cellulase, amino acids) or cheats
    #A Life Strategy role is also chosen as copiotroph or oligotroph. This controls the max growth rate and transporter affinity of species
    
    BQtraits <- 
      runif(n = 20, min = 0, max = 1)
    LStraits <- 
      sample(c(0,1),20,replace=T)
    
    # empty list of speices traits
    S <- 
      data.frame('Species' = c('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11', 'S12',
                                  'S13', 'S14', 'S15', 'S16', 'S17', 'S18', 'S19', 'S20'), 
                 'BQ' = rep(NA, 20), 
                 'LS' = rep(NA, 20), 
                 'Cellulolytic' = rep(F, 20), 
                 'Prototroph' = rep(F, 20))
    # provide the 20 species in S with traits 
    for(i in 1:20){
      Si <- BQtraits[i]
      if(Si < 0.25){
        S$BQ[i] <- 'Yy'
      }else if(Si > 0.25 & Si < 0.5){
        S$BQ[i] <- 'Yn'
      }else if(Si > 0.5 & Si < 0.75){
        S$BQ[i] <- 'Ny'
      }else{
        S$BQ[i] <- 'Nn'
      }
    }
    for(i in 1:20){
      Si <- BQtraits[i]
      if(Si < 0.25){
        S$Cellulolytic[i] <- T
      }else if(Si > 0.5 & Si < 0.75){
        S$Cellulolytic[i] <- T
      }
    }
    for(i in 1:20){
      Si <- BQtraits[i]
      if(Si < 0.5){
        S$Prototroph[i] <- T
      }
    }
    for(i in 1:20){
      Si <- LStraits[i]
      if(Si == 0){
        S$LS[i] <- 'Copio'
      }else{
        S$LS[i] <- 'Oligo'
      }
    }
    
    
    #Set up data frames for Storing proportionality informations, first line shows initial conditions
    #Phototrophs
    Cell.Prototrophs.P$CP.Copios.Prop[y] <- 
      sum(S$BQ %in% 'Yy' & S$LS %in% 'Copio')/20
    Cell.Prototrophs.P$CP.Oligos.Prop[y] <- 
      sum(S$BQ %in% 'Yy' & S$LS %in% 'Oligo')/20
    Cell.Prototrophs.P$CP.Prop[y] <- 
      Cell.Prototrophs.P$CP.Copios.Prop[y] + Cell.Prototrophs.P$CP.Oligos.Prop[y]
    NonCell.Prototrophs.P$NCP.Copios.Prop[y] <- 
      sum(S$BQ %in% 'Yn' & S$LS %in% 'Copio')/20
    NonCell.Prototrophs.P$NCP.Oligos.Prop[y] <- 
      sum(S$BQ %in% 'Yn' & S$LS %in% 'Oligo')/20
    NonCell.Prototrophs.P$NCP.Prop[y] <- 
      NonCell.Prototrophs.P$NCP.Copios.Prop[y] + NonCell.Prototrophs.P$NCP.Oligos.Prop[y]
    #Auxotrophs
    Cell.Auxotrophs.P$CA.Copios.Prop[y] <- 
      sum(S$BQ %in% 'Ny' & S$LS %in% 'Copio')/20
    Cell.Auxotrophs.P$CA.Oligos.Prop[y] <- 
      sum(S$BQ %in% 'Ny' & S$LS %in% 'Oligo')/20
    Cell.Auxotrophs.P$CA.Prop[y] <- 
      Cell.Auxotrophs.P$CA.Copios.Prop[y] + Cell.Auxotrophs.P$CA.Oligos.Prop[y]
    NonCell.Auxotrophs.P$NCA.Copios.Prop[y] <- 
      sum(S$BQ %in% 'Nn' & S$LS %in% 'Copio')/20
    NonCell.Auxotrophs.P$NCA.Oligos.Prop[y] <- 
      sum(S$BQ %in% 'Nn' & S$LS %in% 'Oligo')/20
    NonCell.Auxotrophs.P$NCA.Prop[y] <- 
      NonCell.Auxotrophs.P$NCA.Copios.Prop[y] + NonCell.Auxotrophs.P$NCA.Oligos.Prop[y]
    
    #Step 2: Create params for ODE ----
    
    ## set up constants ----
    times <- 
      seq(0, 500, by = 10) #Time series (generations)
    pars <- 
      c(alpha = 8, gamma = 0.2, Cumax = 0.02, CKm = 1, Oumax = 0.002, OKm = 0.1, 
        CCP = 0.92, CCA = 0.9, CNP = 0.84, CNA = 0.82, 
        OCP = 0.94, OCA = 0.92, ONP = 0.86, ONA = 0.84)
    ## Alpha = glucose units per cellulose;
    ## gamma = death rate per generation (%);
    ## Cumax = max growth rate;
    ## Ckm = half saturation
    ## constant for substrate uptake;
    ## CellPmain, CellAmain, Pmain, NCAmain = relative maintenance energy for each BQ species
    
    ##Startvalues----
    init <- 
      c(C = 2000, G = 100, A = 1e-03, S1 = 0.001, S2 = 0.001, S3 = 0.001, S4 = 0.001, S5 = 0.001, S6 = 0.001, S7 = 0.001,
        S8 = 0.001, S9 = 0.001, S10 = 0.001, S11 = 0.001, S12 = 0.001, S13 = 0.001, S14 = 0.001, S15 = 0.001,
        S16 = 0.001, S17 = 0.001, S18 = 0.001, S19 = 0.001, S20 = 0.001)
    #Initial inputs:
    #C = 2000 uM cellulose;
    #G = 100 uM glucose; 
    #A = 0.001 nM amino acids; 
    #Si = 0.001 ng biomass
    
    ## Set up the ODEs of first timestep ----
    print("Set up initial ODEs")
    tmp <- 
      SetUpODEs(S,BulkSoil) #This configures the initial ODEs, which are a chaotic mess of characters and numerics. Must be eval(parse(...)) for R to read
    
    if(BulkSoil == T){
      Cellevent <- 
        data.frame('var' = rep("C", 4), 'time' = c(100, 200, 300, 400), 
                   'value' = rep(2000, 4), 'method' = rep("add", 4)) #2000 uM cellulose added every 100 generations
    }else if(BulkSoil == F){
      Cellevent <-
        data.frame('var' = rep("G", 50), 'time' = seq(10, 500, by = 10), 
                   'value' = rep(2000, 50), 'method' = rep("add", 50)) 
    }
    
    
    # Step 3: Create list of random mutation events over time, where species lose capacity to generate public goods cellulase oder amino acids ----
    # Must order based on which species mutates first to last
    print("Set up mutaion events")
    
    # empty frame for when mutations occure
    mutationtimes <- 
      rep(0, 20)
    
   ##We create a series of random mutation events among the 10 species----
    mutationtimes <- 
      mutationevents(mutationtimes, MutationRate) 
    
    mutationtimes <- 
      cbind.data.frame('Species' = c('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10',
                                     'S11', 'S12', 'S13', 'S14', 'S15', 'S16', 'S17',
                                     'S18', 'S19', 'S20'), 
                       mutationtimes)
    
    mutationtimes <- 
      mutationtimes[order(mutationtimes$mutationtimes, decreasing = F),]
    mutationcounts$Mutations[y] <- 
      sum(mutationtimes$mutationtimes > 0)
    print(paste0("At least one mutation occures at timesteps ",(paste0(unique(mutationtimes$mutationtimes),sep=" ", collapse = ""))))
    print(paste0("Sum of mutations ", mutationcounts$Mutations[y]))
    
    #Step 4: Create updated lists of mutated species and run SetUpODEs function on the replacement species.---- 
    #This must be done for each time a mutation occurs, as the updated events must build upon each other as the time series progresses
    print("Initialise mutation profils")
    mutantS <- S[,]
    mutantS1 <- mutationsprofile(mutantS, mutationtimes, 1)
    mutantS2 <- mutationsprofile(mutantS1, mutationtimes, 2)
    mutantS3 <- mutationsprofile(mutantS2, mutationtimes, 3)
    mutantS4 <- mutationsprofile(mutantS3, mutationtimes, 4)
    mutantS5 <- mutationsprofile(mutantS4, mutationtimes, 5)
    mutantS6 <- mutationsprofile(mutantS5, mutationtimes, 6)
    mutantS7 <- mutationsprofile(mutantS6, mutationtimes, 7)
    mutantS8 <- mutationsprofile(mutantS7, mutationtimes, 8)
    mutantS9 <- mutationsprofile(mutantS8, mutationtimes, 9)
    mutantS10 <- mutationsprofile(mutantS9, mutationtimes, 10)
    mutantS11 <- mutationsprofile(mutantS10, mutationtimes, 11)
    mutantS12 <- mutationsprofile(mutantS11, mutationtimes, 12)
    mutantS13 <- mutationsprofile(mutantS12, mutationtimes, 13)
    mutantS14 <- mutationsprofile(mutantS13, mutationtimes, 14)
    mutantS15 <- mutationsprofile(mutantS14, mutationtimes, 15)
    mutantS16 <- mutationsprofile(mutantS15, mutationtimes, 16)
    mutantS17 <- mutationsprofile(mutantS16, mutationtimes, 17)
    mutantS18 <- mutationsprofile(mutantS17, mutationtimes, 18)
    mutantS19 <- mutationsprofile(mutantS18, mutationtimes, 19)
    mutantS20 <- mutationsprofile(mutantS19, mutationtimes, 20)
    
    print("Set up mutant ODEs")
    E1 <- SetUpODEs(mutantS1,BulkSoil)
    E2 <- SetUpODEs(mutantS2,BulkSoil)
    E3 <- SetUpODEs(mutantS3,BulkSoil)
    E4 <- SetUpODEs(mutantS4,BulkSoil)
    E5 <- SetUpODEs(mutantS5,BulkSoil)
    E6 <- SetUpODEs(mutantS6,BulkSoil)
    E7 <- SetUpODEs(mutantS7,BulkSoil)
    E8 <- SetUpODEs(mutantS8,BulkSoil)
    E9 <- SetUpODEs(mutantS9,BulkSoil)
    E10 <- SetUpODEs(mutantS10,BulkSoil)
    E11 <- SetUpODEs(mutantS11,BulkSoil)
    E12 <- SetUpODEs(mutantS12,BulkSoil)
    E13 <- SetUpODEs(mutantS13,BulkSoil)
    E14 <- SetUpODEs(mutantS14,BulkSoil)
    E15 <- SetUpODEs(mutantS15,BulkSoil)
    E16 <- SetUpODEs(mutantS16,BulkSoil)
    E17 <- SetUpODEs(mutantS17,BulkSoil)
    E18 <- SetUpODEs(mutantS18,BulkSoil)
    E19 <- SetUpODEs(mutantS19,BulkSoil)
    E20 <- SetUpODEs(mutantS20,BulkSoil)
    
    #Step 5: Run the ODEs
    print("Running ODEs")
    
    # change ODEs for each timestep of those were mutations occured
    CelluloseDegSim <- function(t, state, pars){ #Backbone ODEs, with starting species , after each time step the ODEs are switched with mutation ODEs are inser
        with(as.list(c(state, pars)), {
          d_C <- eval(parse(text=tmp[[1]]))
          d_G <- eval(parse(text=tmp[[2]]))
          d_A <- eval(parse(text=tmp[[3]]))
          d_S1 <- eval(parse(text=tmp$d_Si[1]))
          d_S2 <- eval(parse(text=tmp$d_Si[2]))
          d_S3 <- eval(parse(text=tmp$d_Si[3]))
          d_S4 <- eval(parse(text=tmp$d_Si[4]))
          d_S5 <- eval(parse(text=tmp$d_Si[5]))
          d_S6 <- eval(parse(text=tmp$d_Si[6]))
          d_S7 <- eval(parse(text=tmp$d_Si[7]))
          d_S8 <- eval(parse(text=tmp$d_Si[8]))
          d_S9 <- eval(parse(text=tmp$d_Si[9]))
          d_S10 <- eval(parse(text=tmp$d_Si[10]))
          d_S11 <- eval(parse(text=tmp$d_Si[11]))
          d_S12 <- eval(parse(text=tmp$d_Si[12]))
          d_S13 <- eval(parse(text=tmp$d_Si[13]))
          d_S14 <- eval(parse(text=tmp$d_Si[14]))
          d_S15 <- eval(parse(text=tmp$d_Si[15]))
          d_S16 <- eval(parse(text=tmp$d_Si[16]))
          d_S17 <- eval(parse(text=tmp$d_Si[17]))
          d_S18 <- eval(parse(text=tmp$d_Si[18]))
          d_S19 <- eval(parse(text=tmp$d_Si[19]))
          d_S20 <- eval(parse(text=tmp$d_Si[20]))
          if(t >= mutationtimes$mutationtimes[1]){
            d_C <- eval(parse(text=E1[[1]]))
            d_G <- eval(parse(text=E1[[2]]))
            d_A <- eval(parse(text=E1[[3]]))
            d_S1 <- eval(parse(text=E1$d_Si[1]))
            d_S2 <- eval(parse(text=E1$d_Si[2]))
            d_S3 <- eval(parse(text=E1$d_Si[3]))
            d_S4 <- eval(parse(text=E1$d_Si[4]))
            d_S5 <- eval(parse(text=E1$d_Si[5]))
            d_S6 <- eval(parse(text=E1$d_Si[6]))
            d_S7 <- eval(parse(text=E1$d_Si[7]))
            d_S8 <- eval(parse(text=E1$d_Si[8]))
            d_S9 <- eval(parse(text=E1$d_Si[9]))
            d_S10 <- eval(parse(text=E1$d_Si[10]))
            d_S11 <- eval(parse(text=E1$d_Si[11]))
            d_S12 <- eval(parse(text=E1$d_Si[12]))
            d_S13 <- eval(parse(text=E1$d_Si[13]))
            d_S14 <- eval(parse(text=E1$d_Si[14]))
            d_S15 <- eval(parse(text=E1$d_Si[15]))
            d_S16 <- eval(parse(text=E1$d_Si[16]))
            d_S17 <- eval(parse(text=E1$d_Si[17]))
            d_S18 <- eval(parse(text=E1$d_Si[18]))
            d_S19 <- eval(parse(text=E1$d_Si[19]))
            d_S20 <- eval(parse(text=E1$d_Si[20]))
          }
          if(t >= mutationtimes$mutationtimes[2]){
            d_C <- eval(parse(text=E2[[1]]))
            d_G <- eval(parse(text=E2[[2]]))
            d_A <- eval(parse(text=E2[[3]]))
            d_S1 <- eval(parse(text=E2$d_Si[1]))
            d_S2 <- eval(parse(text=E2$d_Si[2]))
            d_S3 <- eval(parse(text=E2$d_Si[3]))
            d_S4 <- eval(parse(text=E2$d_Si[4]))
            d_S5 <- eval(parse(text=E2$d_Si[5]))
            d_S6 <- eval(parse(text=E2$d_Si[6]))
            d_S7 <- eval(parse(text=E2$d_Si[7]))
            d_S8 <- eval(parse(text=E2$d_Si[8]))
            d_S9 <- eval(parse(text=E2$d_Si[9]))
            d_S10 <- eval(parse(text=E2$d_Si[10]))
            d_S11 <- eval(parse(text=E2$d_Si[11]))
            d_S12 <- eval(parse(text=E2$d_Si[12]))
            d_S13 <- eval(parse(text=E2$d_Si[13]))
            d_S14 <- eval(parse(text=E2$d_Si[14]))
            d_S15 <- eval(parse(text=E2$d_Si[15]))
            d_S16 <- eval(parse(text=E2$d_Si[16]))
            d_S17 <- eval(parse(text=E2$d_Si[17]))
            d_S18 <- eval(parse(text=E2$d_Si[18]))
            d_S19 <- eval(parse(text=E2$d_Si[19]))
            d_S20 <- eval(parse(text=E2$d_Si[20]))
          }
          if(t >= mutationtimes$mutationtimes[3]){
            d_C <- eval(parse(text=E3[[1]]))
            d_G <- eval(parse(text=E3[[2]]))
            d_A <- eval(parse(text=E3[[3]]))
            d_S1 <- eval(parse(text=E3$d_Si[1]))
            d_S2 <- eval(parse(text=E3$d_Si[2]))
            d_S3 <- eval(parse(text=E3$d_Si[3]))
            d_S4 <- eval(parse(text=E3$d_Si[4]))
            d_S5 <- eval(parse(text=E3$d_Si[5]))
            d_S6 <- eval(parse(text=E3$d_Si[6]))
            d_S7 <- eval(parse(text=E3$d_Si[7]))
            d_S8 <- eval(parse(text=E3$d_Si[8]))
            d_S9 <- eval(parse(text=E3$d_Si[9]))
            d_S10 <- eval(parse(text=E3$d_Si[10]))
            d_S11 <- eval(parse(text=E3$d_Si[11]))
            d_S12 <- eval(parse(text=E3$d_Si[12]))
            d_S13 <- eval(parse(text=E3$d_Si[13]))
            d_S14 <- eval(parse(text=E3$d_Si[14]))
            d_S15 <- eval(parse(text=E3$d_Si[15]))
            d_S16 <- eval(parse(text=E3$d_Si[16]))
            d_S17 <- eval(parse(text=E3$d_Si[17]))
            d_S18 <- eval(parse(text=E3$d_Si[18]))
            d_S19 <- eval(parse(text=E3$d_Si[19]))
            d_S20 <- eval(parse(text=E3$d_Si[20]))
          }
          if(t >= mutationtimes$mutationtimes[4]){
            d_C <- eval(parse(text=E4[[1]]))
            d_G <- eval(parse(text=E4[[2]]))
            d_A <- eval(parse(text=E4[[3]]))
            d_S1 <- eval(parse(text=E4$d_Si[1]))
            d_S2 <- eval(parse(text=E4$d_Si[2]))
            d_S3 <- eval(parse(text=E4$d_Si[3]))
            d_S4 <- eval(parse(text=E4$d_Si[4]))
            d_S5 <- eval(parse(text=E4$d_Si[5]))
            d_S6 <- eval(parse(text=E4$d_Si[6]))
            d_S7 <- eval(parse(text=E4$d_Si[7]))
            d_S8 <- eval(parse(text=E4$d_Si[8]))
            d_S9 <- eval(parse(text=E4$d_Si[9]))
            d_S10 <- eval(parse(text=E4$d_Si[10]))
            d_S11 <- eval(parse(text=E4$d_Si[11]))
            d_S12 <- eval(parse(text=E4$d_Si[12]))
            d_S13 <- eval(parse(text=E4$d_Si[13]))
            d_S14 <- eval(parse(text=E4$d_Si[14]))
            d_S15 <- eval(parse(text=E4$d_Si[15]))
            d_S16 <- eval(parse(text=E4$d_Si[16]))
            d_S17 <- eval(parse(text=E4$d_Si[17]))
            d_S18 <- eval(parse(text=E4$d_Si[18]))
            d_S19 <- eval(parse(text=E4$d_Si[19]))
            d_S20 <- eval(parse(text=E4$d_Si[20]))
          }
          if(t >= mutationtimes$mutationtimes[5]){
            d_C <- eval(parse(text=E5[[1]]))
            d_G <- eval(parse(text=E5[[2]]))
            d_A <- eval(parse(text=E5[[3]]))
            d_S1 <- eval(parse(text=E5$d_Si[1]))
            d_S2 <- eval(parse(text=E5$d_Si[2]))
            d_S3 <- eval(parse(text=E5$d_Si[3]))
            d_S4 <- eval(parse(text=E5$d_Si[4]))
            d_S5 <- eval(parse(text=E5$d_Si[5]))
            d_S6 <- eval(parse(text=E5$d_Si[6]))
            d_S7 <- eval(parse(text=E5$d_Si[7]))
            d_S8 <- eval(parse(text=E5$d_Si[8]))
            d_S9 <- eval(parse(text=E5$d_Si[9]))
            d_S10 <- eval(parse(text=E5$d_Si[10]))
            d_S11 <- eval(parse(text=E5$d_Si[11]))
            d_S12 <- eval(parse(text=E5$d_Si[12]))
            d_S13 <- eval(parse(text=E5$d_Si[13]))
            d_S14 <- eval(parse(text=E5$d_Si[14]))
            d_S15 <- eval(parse(text=E5$d_Si[15]))
            d_S16 <- eval(parse(text=E5$d_Si[16]))
            d_S17 <- eval(parse(text=E5$d_Si[17]))
            d_S18 <- eval(parse(text=E5$d_Si[18]))
            d_S19 <- eval(parse(text=E5$d_Si[19]))
            d_S20 <- eval(parse(text=E5$d_Si[20]))
          }
          if(t >= mutationtimes$mutationtimes[6]){
            d_C <- eval(parse(text=E6[[1]]))
            d_G <- eval(parse(text=E6[[2]]))
            d_A <- eval(parse(text=E6[[3]]))
            d_S1 <- eval(parse(text=E6$d_Si[1]))
            d_S2 <- eval(parse(text=E6$d_Si[2]))
            d_S3 <- eval(parse(text=E6$d_Si[3]))
            d_S4 <- eval(parse(text=E6$d_Si[4]))
            d_S5 <- eval(parse(text=E6$d_Si[5]))
            d_S6 <- eval(parse(text=E6$d_Si[6]))
            d_S7 <- eval(parse(text=E6$d_Si[7]))
            d_S8 <- eval(parse(text=E6$d_Si[8]))
            d_S9 <- eval(parse(text=E6$d_Si[9]))
            d_S10 <- eval(parse(text=E6$d_Si[10]))
            d_S11 <- eval(parse(text=E6$d_Si[11]))
            d_S12 <- eval(parse(text=E6$d_Si[12]))
            d_S13 <- eval(parse(text=E6$d_Si[13]))
            d_S14 <- eval(parse(text=E6$d_Si[14]))
            d_S15 <- eval(parse(text=E6$d_Si[15]))
            d_S16 <- eval(parse(text=E6$d_Si[16]))
            d_S17 <- eval(parse(text=E6$d_Si[17]))
            d_S18 <- eval(parse(text=E6$d_Si[18]))
            d_S19 <- eval(parse(text=E6$d_Si[19]))
            d_S20 <- eval(parse(text=E6$d_Si[20]))
          }
          if(t >= mutationtimes$mutationtimes[7]){
            d_C <- eval(parse(text=E7[[1]]))
            d_G <- eval(parse(text=E7[[2]]))
            d_A <- eval(parse(text=E7[[3]]))
            d_S1 <- eval(parse(text=E7$d_Si[1]))
            d_S2 <- eval(parse(text=E7$d_Si[2]))
            d_S3 <- eval(parse(text=E7$d_Si[3]))
            d_S4 <- eval(parse(text=E7$d_Si[4]))
            d_S5 <- eval(parse(text=E7$d_Si[5]))
            d_S6 <- eval(parse(text=E7$d_Si[6]))
            d_S7 <- eval(parse(text=E7$d_Si[7]))
            d_S8 <- eval(parse(text=E7$d_Si[8]))
            d_S9 <- eval(parse(text=E7$d_Si[9]))
            d_S10 <- eval(parse(text=E7$d_Si[10]))
            d_S11 <- eval(parse(text=E7$d_Si[11]))
            d_S12 <- eval(parse(text=E7$d_Si[12]))
            d_S13 <- eval(parse(text=E7$d_Si[13]))
            d_S14 <- eval(parse(text=E7$d_Si[14]))
            d_S15 <- eval(parse(text=E7$d_Si[15]))
            d_S16 <- eval(parse(text=E7$d_Si[16]))
            d_S17 <- eval(parse(text=E7$d_Si[17]))
            d_S18 <- eval(parse(text=E7$d_Si[18]))
            d_S19 <- eval(parse(text=E7$d_Si[19]))
            d_S20 <- eval(parse(text=E7$d_Si[20]))
          }
          if(t >= mutationtimes$mutationtimes[8]){
            d_C <- eval(parse(text=E8[[1]]))
            d_G <- eval(parse(text=E8[[2]]))
            d_A <- eval(parse(text=E8[[3]]))
            d_S1 <- eval(parse(text=E8$d_Si[1]))
            d_S2 <- eval(parse(text=E8$d_Si[2]))
            d_S3 <- eval(parse(text=E8$d_Si[3]))
            d_S4 <- eval(parse(text=E8$d_Si[4]))
            d_S5 <- eval(parse(text=E8$d_Si[5]))
            d_S6 <- eval(parse(text=E8$d_Si[6]))
            d_S7 <- eval(parse(text=E8$d_Si[7]))
            d_S8 <- eval(parse(text=E8$d_Si[8]))
            d_S9 <- eval(parse(text=E8$d_Si[9]))
            d_S10 <- eval(parse(text=E8$d_Si[10]))
            d_S11 <- eval(parse(text=E8$d_Si[11]))
            d_S12 <- eval(parse(text=E8$d_Si[12]))
            d_S13 <- eval(parse(text=E8$d_Si[13]))
            d_S14 <- eval(parse(text=E8$d_Si[14]))
            d_S15 <- eval(parse(text=E8$d_Si[15]))
            d_S16 <- eval(parse(text=E8$d_Si[16]))
            d_S17 <- eval(parse(text=E8$d_Si[17]))
            d_S18 <- eval(parse(text=E8$d_Si[18]))
            d_S19 <- eval(parse(text=E8$d_Si[19]))
            d_S20 <- eval(parse(text=E8$d_Si[20]))
          }
          if(t >= mutationtimes$mutationtimes[9]){
            d_C <- eval(parse(text=E9[[1]]))
            d_G <- eval(parse(text=E9[[2]]))
            d_A <- eval(parse(text=E9[[3]]))
            d_S1 <- eval(parse(text=E9$d_Si[1]))
            d_S2 <- eval(parse(text=E9$d_Si[2]))
            d_S3 <- eval(parse(text=E9$d_Si[3]))
            d_S4 <- eval(parse(text=E9$d_Si[4]))
            d_S5 <- eval(parse(text=E9$d_Si[5]))
            d_S6 <- eval(parse(text=E9$d_Si[6]))
            d_S7 <- eval(parse(text=E9$d_Si[7]))
            d_S8 <- eval(parse(text=E9$d_Si[8]))
            d_S9 <- eval(parse(text=E9$d_Si[9]))
            d_S10 <- eval(parse(text=E9$d_Si[10]))
            d_S11 <- eval(parse(text=E9$d_Si[11]))
            d_S12 <- eval(parse(text=E9$d_Si[12]))
            d_S13 <- eval(parse(text=E9$d_Si[13]))
            d_S14 <- eval(parse(text=E9$d_Si[14]))
            d_S15 <- eval(parse(text=E9$d_Si[15]))
            d_S16 <- eval(parse(text=E9$d_Si[16]))
            d_S17 <- eval(parse(text=E9$d_Si[17]))
            d_S18 <- eval(parse(text=E9$d_Si[18]))
            d_S19 <- eval(parse(text=E9$d_Si[19]))
            d_S20 <- eval(parse(text=E9$d_Si[20]))
          }
          if(t >= mutationtimes$mutationtimes[10]){
            d_C <- eval(parse(text=E10[[1]]))
            d_G <- eval(parse(text=E10[[2]]))
            d_A <- eval(parse(text=E10[[3]]))
            d_S1 <- eval(parse(text=E10$d_Si[1]))
            d_S2 <- eval(parse(text=E10$d_Si[2]))
            d_S3 <- eval(parse(text=E10$d_Si[3]))
            d_S4 <- eval(parse(text=E10$d_Si[4]))
            d_S5 <- eval(parse(text=E10$d_Si[5]))
            d_S6 <- eval(parse(text=E10$d_Si[6]))
            d_S7 <- eval(parse(text=E10$d_Si[7]))
            d_S8 <- eval(parse(text=E10$d_Si[8]))
            d_S9 <- eval(parse(text=E10$d_Si[9]))
            d_S10 <- eval(parse(text=E10$d_Si[10]))
            d_S11 <- eval(parse(text=E10$d_Si[11]))
            d_S12 <- eval(parse(text=E10$d_Si[12]))
            d_S13 <- eval(parse(text=E10$d_Si[13]))
            d_S14 <- eval(parse(text=E10$d_Si[14]))
            d_S15 <- eval(parse(text=E10$d_Si[15]))
            d_S16 <- eval(parse(text=E10$d_Si[16]))
            d_S17 <- eval(parse(text=E10$d_Si[17]))
            d_S18 <- eval(parse(text=E10$d_Si[18]))
            d_S19 <- eval(parse(text=E10$d_Si[19]))
            d_S20 <- eval(parse(text=E10$d_Si[20]))
          }
          if(t >= mutationtimes$mutationtimes[11]){
            d_C <- eval(parse(text=E11[[1]]))
            d_G <- eval(parse(text=E11[[2]]))
            d_A <- eval(parse(text=E11[[3]]))
            d_S1 <- eval(parse(text=E11$d_Si[1]))
            d_S2 <- eval(parse(text=E11$d_Si[2]))
            d_S3 <- eval(parse(text=E11$d_Si[3]))
            d_S4 <- eval(parse(text=E11$d_Si[4]))
            d_S5 <- eval(parse(text=E11$d_Si[5]))
            d_S6 <- eval(parse(text=E11$d_Si[6]))
            d_S7 <- eval(parse(text=E11$d_Si[7]))
            d_S8 <- eval(parse(text=E11$d_Si[8]))
            d_S9 <- eval(parse(text=E11$d_Si[9]))
            d_S10 <- eval(parse(text=E11$d_Si[10]))
            d_S11 <- eval(parse(text=E11$d_Si[11]))
            d_S12 <- eval(parse(text=E11$d_Si[12]))
            d_S13 <- eval(parse(text=E11$d_Si[13]))
            d_S14 <- eval(parse(text=E11$d_Si[14]))
            d_S15 <- eval(parse(text=E11$d_Si[15]))
            d_S16 <- eval(parse(text=E11$d_Si[16]))
            d_S17 <- eval(parse(text=E11$d_Si[17]))
            d_S18 <- eval(parse(text=E11$d_Si[18]))
            d_S19 <- eval(parse(text=E11$d_Si[19]))
            d_S20 <- eval(parse(text=E11$d_Si[20]))
          }
          if(t >= mutationtimes$mutationtimes[12]){
            d_C <- eval(parse(text=E12[[1]]))
            d_G <- eval(parse(text=E12[[2]]))
            d_A <- eval(parse(text=E12[[3]]))
            d_S1 <- eval(parse(text=E12$d_Si[1]))
            d_S2 <- eval(parse(text=E12$d_Si[2]))
            d_S3 <- eval(parse(text=E12$d_Si[3]))
            d_S4 <- eval(parse(text=E12$d_Si[4]))
            d_S5 <- eval(parse(text=E12$d_Si[5]))
            d_S6 <- eval(parse(text=E12$d_Si[6]))
            d_S7 <- eval(parse(text=E12$d_Si[7]))
            d_S8 <- eval(parse(text=E12$d_Si[8]))
            d_S9 <- eval(parse(text=E12$d_Si[9]))
            d_S10 <- eval(parse(text=E12$d_Si[10]))
            d_S11 <- eval(parse(text=E12$d_Si[11]))
            d_S12 <- eval(parse(text=E12$d_Si[12]))
            d_S13 <- eval(parse(text=E12$d_Si[13]))
            d_S14 <- eval(parse(text=E12$d_Si[14]))
            d_S15 <- eval(parse(text=E12$d_Si[15]))
            d_S16 <- eval(parse(text=E12$d_Si[16]))
            d_S17 <- eval(parse(text=E12$d_Si[17]))
            d_S18 <- eval(parse(text=E12$d_Si[18]))
            d_S19 <- eval(parse(text=E12$d_Si[19]))
            d_S20 <- eval(parse(text=E12$d_Si[20]))
          }
          if(t >= mutationtimes$mutationtimes[13]){
            d_C <- eval(parse(text=E13[[1]]))
            d_G <- eval(parse(text=E13[[2]]))
            d_A <- eval(parse(text=E13[[3]]))
            d_S1 <- eval(parse(text=E13$d_Si[1]))
            d_S2 <- eval(parse(text=E13$d_Si[2]))
            d_S3 <- eval(parse(text=E13$d_Si[3]))
            d_S4 <- eval(parse(text=E13$d_Si[4]))
            d_S5 <- eval(parse(text=E13$d_Si[5]))
            d_S6 <- eval(parse(text=E13$d_Si[6]))
            d_S7 <- eval(parse(text=E13$d_Si[7]))
            d_S8 <- eval(parse(text=E13$d_Si[8]))
            d_S9 <- eval(parse(text=E13$d_Si[9]))
            d_S10 <- eval(parse(text=E13$d_Si[10]))
            d_S11 <- eval(parse(text=E13$d_Si[11]))
            d_S12 <- eval(parse(text=E13$d_Si[12]))
            d_S13 <- eval(parse(text=E13$d_Si[13]))
            d_S14 <- eval(parse(text=E13$d_Si[14]))
            d_S15 <- eval(parse(text=E13$d_Si[15]))
            d_S16 <- eval(parse(text=E13$d_Si[16]))
            d_S17 <- eval(parse(text=E13$d_Si[17]))
            d_S18 <- eval(parse(text=E13$d_Si[18]))
            d_S19 <- eval(parse(text=E13$d_Si[19]))
            d_S20 <- eval(parse(text=E13$d_Si[20]))
          }
          if(t >= mutationtimes$mutationtimes[14]){
            d_C <- eval(parse(text=E14[[1]]))
            d_G <- eval(parse(text=E14[[2]]))
            d_A <- eval(parse(text=E14[[3]]))
            d_S1 <- eval(parse(text=E14$d_Si[1]))
            d_S2 <- eval(parse(text=E14$d_Si[2]))
            d_S3 <- eval(parse(text=E14$d_Si[3]))
            d_S4 <- eval(parse(text=E14$d_Si[4]))
            d_S5 <- eval(parse(text=E14$d_Si[5]))
            d_S6 <- eval(parse(text=E14$d_Si[6]))
            d_S7 <- eval(parse(text=E14$d_Si[7]))
            d_S8 <- eval(parse(text=E14$d_Si[8]))
            d_S9 <- eval(parse(text=E14$d_Si[9]))
            d_S10 <- eval(parse(text=E14$d_Si[10]))
            d_S11 <- eval(parse(text=E14$d_Si[11]))
            d_S12 <- eval(parse(text=E14$d_Si[12]))
            d_S13 <- eval(parse(text=E14$d_Si[13]))
            d_S14 <- eval(parse(text=E14$d_Si[14]))
            d_S15 <- eval(parse(text=E14$d_Si[15]))
            d_S16 <- eval(parse(text=E14$d_Si[16]))
            d_S17 <- eval(parse(text=E14$d_Si[17]))
            d_S18 <- eval(parse(text=E14$d_Si[18]))
            d_S19 <- eval(parse(text=E14$d_Si[19]))
            d_S20 <- eval(parse(text=E14$d_Si[20]))
          }
          if(t >= mutationtimes$mutationtimes[15]){
            d_C <- eval(parse(text=E15[[1]]))
            d_G <- eval(parse(text=E15[[2]]))
            d_A <- eval(parse(text=E15[[3]]))
            d_S1 <- eval(parse(text=E15$d_Si[1]))
            d_S2 <- eval(parse(text=E15$d_Si[2]))
            d_S3 <- eval(parse(text=E15$d_Si[3]))
            d_S4 <- eval(parse(text=E15$d_Si[4]))
            d_S5 <- eval(parse(text=E15$d_Si[5]))
            d_S6 <- eval(parse(text=E15$d_Si[6]))
            d_S7 <- eval(parse(text=E15$d_Si[7]))
            d_S8 <- eval(parse(text=E15$d_Si[8]))
            d_S9 <- eval(parse(text=E15$d_Si[9]))
            d_S10 <- eval(parse(text=E15$d_Si[10]))
            d_S11 <- eval(parse(text=E15$d_Si[11]))
            d_S12 <- eval(parse(text=E15$d_Si[12]))
            d_S13 <- eval(parse(text=E15$d_Si[13]))
            d_S14 <- eval(parse(text=E15$d_Si[14]))
            d_S15 <- eval(parse(text=E15$d_Si[15]))
            d_S16 <- eval(parse(text=E15$d_Si[16]))
            d_S17 <- eval(parse(text=E15$d_Si[17]))
            d_S18 <- eval(parse(text=E15$d_Si[18]))
            d_S19 <- eval(parse(text=E15$d_Si[19]))
            d_S20 <- eval(parse(text=E15$d_Si[20]))
          }
          if(t >= mutationtimes$mutationtimes[16]){
            d_C <- eval(parse(text=E16[[1]]))
            d_G <- eval(parse(text=E16[[2]]))
            d_A <- eval(parse(text=E16[[3]]))
            d_S1 <- eval(parse(text=E16$d_Si[1]))
            d_S2 <- eval(parse(text=E16$d_Si[2]))
            d_S3 <- eval(parse(text=E16$d_Si[3]))
            d_S4 <- eval(parse(text=E16$d_Si[4]))
            d_S5 <- eval(parse(text=E16$d_Si[5]))
            d_S6 <- eval(parse(text=E16$d_Si[6]))
            d_S7 <- eval(parse(text=E16$d_Si[7]))
            d_S8 <- eval(parse(text=E16$d_Si[8]))
            d_S9 <- eval(parse(text=E16$d_Si[9]))
            d_S10 <- eval(parse(text=E16$d_Si[10]))
            d_S11 <- eval(parse(text=E16$d_Si[11]))
            d_S12 <- eval(parse(text=E16$d_Si[12]))
            d_S13 <- eval(parse(text=E16$d_Si[13]))
            d_S14 <- eval(parse(text=E16$d_Si[14]))
            d_S15 <- eval(parse(text=E16$d_Si[15]))
            d_S16 <- eval(parse(text=E16$d_Si[16]))
            d_S17 <- eval(parse(text=E16$d_Si[17]))
            d_S18 <- eval(parse(text=E16$d_Si[18]))
            d_S19 <- eval(parse(text=E16$d_Si[19]))
            d_S20 <- eval(parse(text=E16$d_Si[20]))
          }
          if(t >= mutationtimes$mutationtimes[17]){
            d_C <- eval(parse(text=E17[[1]]))
            d_G <- eval(parse(text=E17[[2]]))
            d_A <- eval(parse(text=E17[[3]]))
            d_S1 <- eval(parse(text=E17$d_Si[1]))
            d_S2 <- eval(parse(text=E17$d_Si[2]))
            d_S3 <- eval(parse(text=E17$d_Si[3]))
            d_S4 <- eval(parse(text=E17$d_Si[4]))
            d_S5 <- eval(parse(text=E17$d_Si[5]))
            d_S6 <- eval(parse(text=E17$d_Si[6]))
            d_S7 <- eval(parse(text=E17$d_Si[7]))
            d_S8 <- eval(parse(text=E17$d_Si[8]))
            d_S9 <- eval(parse(text=E17$d_Si[9]))
            d_S10 <- eval(parse(text=E17$d_Si[10]))
            d_S11 <- eval(parse(text=E17$d_Si[11]))
            d_S12 <- eval(parse(text=E17$d_Si[12]))
            d_S13 <- eval(parse(text=E17$d_Si[13]))
            d_S14 <- eval(parse(text=E17$d_Si[14]))
            d_S15 <- eval(parse(text=E17$d_Si[15]))
            d_S16 <- eval(parse(text=E17$d_Si[16]))
            d_S17 <- eval(parse(text=E17$d_Si[17]))
            d_S18 <- eval(parse(text=E17$d_Si[18]))
            d_S19 <- eval(parse(text=E17$d_Si[19]))
            d_S20 <- eval(parse(text=E17$d_Si[20]))
          }
          if(t >= mutationtimes$mutationtimes[18]){
            d_C <- eval(parse(text=E18[[1]]))
            d_G <- eval(parse(text=E18[[2]]))
            d_A <- eval(parse(text=E18[[3]]))
            d_S1 <- eval(parse(text=E18$d_Si[1]))
            d_S2 <- eval(parse(text=E18$d_Si[2]))
            d_S3 <- eval(parse(text=E18$d_Si[3]))
            d_S4 <- eval(parse(text=E18$d_Si[4]))
            d_S5 <- eval(parse(text=E18$d_Si[5]))
            d_S6 <- eval(parse(text=E18$d_Si[6]))
            d_S7 <- eval(parse(text=E18$d_Si[7]))
            d_S8 <- eval(parse(text=E18$d_Si[8]))
            d_S9 <- eval(parse(text=E18$d_Si[9]))
            d_S10 <- eval(parse(text=E18$d_Si[10]))
            d_S11 <- eval(parse(text=E18$d_Si[11]))
            d_S12 <- eval(parse(text=E18$d_Si[12]))
            d_S13 <- eval(parse(text=E18$d_Si[13]))
            d_S14 <- eval(parse(text=E18$d_Si[14]))
            d_S15 <- eval(parse(text=E18$d_Si[15]))
            d_S16 <- eval(parse(text=E18$d_Si[16]))
            d_S17 <- eval(parse(text=E18$d_Si[17]))
            d_S18 <- eval(parse(text=E18$d_Si[18]))
            d_S19 <- eval(parse(text=E18$d_Si[19]))
            d_S20 <- eval(parse(text=E18$d_Si[20]))
          }
          if(t >= mutationtimes$mutationtimes[19]){
            d_C <- eval(parse(text=E19[[1]]))
            d_G <- eval(parse(text=E19[[2]]))
            d_A <- eval(parse(text=E19[[3]]))
            d_S1 <- eval(parse(text=E19$d_Si[1]))
            d_S2 <- eval(parse(text=E19$d_Si[2]))
            d_S3 <- eval(parse(text=E19$d_Si[3]))
            d_S4 <- eval(parse(text=E19$d_Si[4]))
            d_S5 <- eval(parse(text=E19$d_Si[5]))
            d_S6 <- eval(parse(text=E19$d_Si[6]))
            d_S7 <- eval(parse(text=E19$d_Si[7]))
            d_S8 <- eval(parse(text=E19$d_Si[8]))
            d_S9 <- eval(parse(text=E19$d_Si[9]))
            d_S10 <- eval(parse(text=E19$d_Si[10]))
            d_S11 <- eval(parse(text=E19$d_Si[11]))
            d_S12 <- eval(parse(text=E19$d_Si[12]))
            d_S13 <- eval(parse(text=E19$d_Si[13]))
            d_S14 <- eval(parse(text=E19$d_Si[14]))
            d_S15 <- eval(parse(text=E19$d_Si[15]))
            d_S16 <- eval(parse(text=E19$d_Si[16]))
            d_S17 <- eval(parse(text=E19$d_Si[17]))
            d_S18 <- eval(parse(text=E19$d_Si[18]))
            d_S19 <- eval(parse(text=E19$d_Si[19]))
            d_S20 <- eval(parse(text=E19$d_Si[20]))
          }
          if(t >= mutationtimes$mutationtimes[20]){
            d_C <- eval(parse(text=E20[[1]]))
            d_G <- eval(parse(text=E20[[2]]))
            d_A <- eval(parse(text=E20[[3]]))
            d_S1 <- eval(parse(text=E20$d_Si[1]))
            d_S2 <- eval(parse(text=E20$d_Si[2]))
            d_S3 <- eval(parse(text=E20$d_Si[3]))
            d_S4 <- eval(parse(text=E20$d_Si[4]))
            d_S5 <- eval(parse(text=E20$d_Si[5]))
            d_S6 <- eval(parse(text=E20$d_Si[6]))
            d_S7 <- eval(parse(text=E20$d_Si[7]))
            d_S8 <- eval(parse(text=E20$d_Si[8]))
            d_S9 <- eval(parse(text=E20$d_Si[9]))
            d_S10 <- eval(parse(text=E20$d_Si[10]))
            d_S11 <- eval(parse(text=E20$d_Si[11]))
            d_S12 <- eval(parse(text=E20$d_Si[12]))
            d_S13 <- eval(parse(text=E20$d_Si[13]))
            d_S14 <- eval(parse(text=E20$d_Si[14]))
            d_S15 <- eval(parse(text=E20$d_Si[15]))
            d_S16 <- eval(parse(text=E20$d_Si[16]))
            d_S17 <- eval(parse(text=E20$d_Si[17]))
            d_S18 <- eval(parse(text=E20$d_Si[18]))
            d_S19 <- eval(parse(text=E20$d_Si[19]))
            d_S20 <- eval(parse(text=E20$d_Si[20]))
          }
          return(list(c(C = d_C, G = d_G, A = d_A, S1 = d_S1, S2 = d_S2, S3 = d_S3, S4 = d_S4, S5 = d_S5, S6 = d_S6, S7 = d_S7, 
                        S8 = d_S8, S9 = d_S9, S10 = d_S10, S11 = d_S11, S12 = d_S12, S13 = d_S13, S14 = d_S14, S15 = d_S15,
                        S16 = d_S16, S17 = d_S17, S18 = d_S18, S19 = d_S19, S20 = d_S20)))
        })
      }

    myres <- 
      ode(init, times, CelluloseDegSim, pars, events = list(data = Cellevent)) ##Executing the ODEs ----
    myres[is.na(myres)] <- 0 #Where species die, NA is the result -- must convert to 0s
    
    #Calculating the cumulative biomass of each species in each simulation
    
    CBSum <- 
      (sum(sum(myres[,5]), sum(myres[,6]), sum(myres[,7]), sum(myres[,8]), 
           sum(myres[,9]), sum(myres[,10]), sum(myres[,11]), sum(myres[,12]),
           sum(myres[,13]), sum(myres[,14]), sum(myres[,15]), sum(myres[,16]),
           sum(myres[,17]), sum(myres[,18]), sum(myres[,19]), sum(myres[,20]),
           sum(myres[,21]) ,sum(myres[,22]), sum(myres[,23]), sum(myres[,24])))
    CBdf$Biomass[y] <- log10(CBSum) #Log10 transform cumulative biomass production 
    
    
    #Keeping track of final species distributions per simulated community
    Scounts[y, ] <- 
      c(myres[,5][nrow(myres)],  myres[,6][nrow(myres)],  
        myres[,7][nrow(myres)],  myres[,8][nrow(myres)], 
        myres[,9][nrow(myres)],  myres[,10][nrow(myres)], 
        myres[,11][nrow(myres)], myres[,12][nrow(myres)], 
        myres[,13][nrow(myres)], myres[,14][nrow(myres)],
        myres[,15][nrow(myres)], myres[,16][nrow(myres)], 
        myres[,17][nrow(myres)], myres[,18][nrow(myres)], 
        myres[,19][nrow(myres)], myres[,20][nrow(myres)], 
        myres[,21][nrow(myres)], myres[,22][nrow(myres)], 
        myres[,23][nrow(myres)], myres[,24][nrow(myres)])
    
    
    #Calculating alpha-diversity for each simulation at the final time point
    
    Shannondf$Alpha.Div[y] <- 
      diversity(Scounts[y,], index = 'shannon')
    Simpsondf$Alpha.Div[y] <- 
      diversity(Scounts[y,], index = 'simpson')
    
    #Calculating the mean cumulative biomass of each species type, and weighting that by the proportion of each species type in each community
    for(x in 1:20){
      Si.CB <- sum(myres[,x+4])
      if(S$BQ[x] %in% 'Yy' & S$LS[x] %in% 'Copio'){
        Cell.Prototrophs.B$CP.Copios.Biom[y] <- 
          Cell.Prototrophs.B$CP.Copios.Biom[y] + (Si.CB*Cell.Prototrophs.P$CP.Copios.Prop[y])
      }else if(S$BQ[x] %in% 'Yy' & S$LS[x] %in% 'Oligo'){
        Cell.Prototrophs.B$CP.Oligos.Biom[y] <- 
          Cell.Prototrophs.B$CP.Oligos.Biom[y] + (Si.CB*Cell.Prototrophs.P$CP.Oligos.Prop[y])
      }else if(S$BQ[x] %in% 'Yn' & S$LS[x] %in% 'Copio'){
        NonCell.Prototrophs.B$NCP.Copios.Biom[y] <- 
          NonCell.Prototrophs.B$NCP.Copios.Biom[y] + (Si.CB*NonCell.Prototrophs.P$NCP.Copios.Prop[y])
      }else if(S$BQ[x] %in% 'Yn' & S$LS[x] %in% 'Oligo'){
        NonCell.Prototrophs.B$NCP.Oligos.Biom[y] <- 
          NonCell.Prototrophs.B$NCP.Oligos.Biom[y] + (Si.CB*NonCell.Prototrophs.P$NCP.Oligos.Prop[y])
      }else if(S$BQ[x] %in% 'Ny' & S$LS[x] %in% 'Copio'){
        Cell.Auxotrophs.B$CA.Copios.Biom[y] <- 
          Cell.Auxotrophs.B$CA.Copios.Biom[y] + (Si.CB*Cell.Auxotrophs.P$CA.Copios.Prop[y])
      }else if(S$BQ[x] %in% 'Ny' & S$LS[x] %in% 'Oligo'){
        Cell.Auxotrophs.B$CA.Oligos.Biom[y] <- 
          Cell.Auxotrophs.B$CA.Oligos.Biom[y] + (Si.CB*Cell.Auxotrophs.P$CA.Oligos.Prop[y])
      }else if(S$BQ[x] %in% 'Nn' & S$LS[x] %in% 'Copio'){
        NonCell.Auxotrophs.B$NCA.Copios.Biom[y] <- NonCell.Auxotrophs.B$NCA.Copios.Biom[y] + (Si.CB*NonCell.Auxotrophs.P$NCA.Copios.Prop[y])
      }else if(S$BQ[x] %in% 'Nn' & S$LS[x] %in% 'Oligo'){
        NonCell.Auxotrophs.B$NCA.Oligos.Biom[y] <- NonCell.Auxotrophs.B$NCA.Oligos.Biom[y] + (Si.CB*NonCell.Auxotrophs.P$NCA.Oligos.Prop[y])
      }
    }
    
    #Simulation ended and information sorted! Now to keep track of everything as it runs, report how long the process takes
    end.time <- Sys.time()
    print('Completed Simulation:')
    print(y)
    print('Processing time:')
    print(end.time - start.time) 
  }
  return(list(CBdf = CBdf, Scounts = Scounts, 
              mutationcounts = mutationcounts,mutationtimes=mutationtimes, 
              Cell.Auxotrophs.B = Cell.Auxotrophs.B, Cell.Auxotrophs.P = Cell.Auxotrophs.P,
              Cell.Prototrophs.B = Cell.Prototrophs.B, Cell.Prototrophs.P=Cell.Prototrophs.P,
              NonCell.Auxotrophs.B = NonCell.Auxotrophs.B, NonCell.Auxotrophs.P=NonCell.Auxotrophs.P,
              NonCell.Prototrophs.B = NonCell.Prototrophs.B, NonCell.Prototrophs.P = NonCell.Prototrophs.P,
              Shannondf = Shannondf, Simpsondf = Simpsondf,
              myres = myres, times = times))
}
