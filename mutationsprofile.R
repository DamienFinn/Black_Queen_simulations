#A function to create new mutation profiles that build upon each other and keep track of mutation over time
#Note that species with two public good functions are randomly mutated to only cellulase or amino acid production
#Species that have only one public good function are mutated to Cheaters


mutationsprofile <- function(mutantS, mutationtimes, eventno){
  if(mutationtimes$mutationtimes[eventno] > 0){
    myref <- mutationtimes$Species[eventno]
    for(i in 1:20){
      mycheck <- mutantS$Species[i]
      checking <- match(mycheck, myref, nomatch = 0)
      if(checking > 0 & mutantS$BQ[i] %in% 'Yy' & mutantS$LS[i] %in% 'Copio'){
        geneloss <- sample(c(0,1),1,replace=T)
        if(geneloss > 0.5){
          mutantS[i,] <- c("M", 'Yn', 'Copio', F, T)
        }else{
          mutantS[i,] <- c("M", 'Ny', 'Copio', T, F)
        }
      }else if(checking > 0 & mutantS$BQ[i] %in% 'Yy' & mutantS$LS[i] %in% 'Oligo'){
        geneloss <- sample(c(0,1),1,replace=T)
        if(geneloss > 0.5){
          mutantS[i,] <- c("M", 'Yn', 'Oligo', F, F)
        }else{
          mutantS[i,] <- c("M", 'Ny', 'Oligo', F, F)
        }
      }else if(checking > 0 & !(mutantS$BQ[i] %in% 'Nn') & mutantS$LS[i] %in% 'Copio'){
        mutantS[i,] <- c("M", 'Nn', 'Copio', F, F)
      }else if(checking > 0 & !(mutantS$BQ[i] %in% 'Nn') & mutantS$LS[i] %in% 'Oligo'){
        mutantS[i,] <- c("M", 'Nn', 'Oligo', F, F)
      }
    }
  }
  return(mutantS)
}  
