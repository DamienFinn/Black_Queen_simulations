#A function to format the changes in cellulose, glucose, amino acids and species biomass in the Cellulose Degradation Simulation
#This must be formatted as characters and then converted back to an expression with eval(parse(...)) when running as ODE

SetUpODEs <- function(S,BulkSoil){
  
  allCellul <- c()
  GlucProd <- c()
  GlucDemand <- c()
  AminoAcidProd <- c()
  AminoAcidDemand <- c()
  d_Sidf <- data.frame('d_Si' = rep(NA, 20))
  
  for(i in 1:20){

      Si <- paste('S', i, sep = '')
      
      if(S$Cellulolytic[i] == T){
        allCellul <- paste(allCellul, '- (', Si, '*19*(C/9.4 + C))', sep = '')
        GlucProd <- paste(GlucProd, '+ alpha*(',Si,'*19*(C/9.4 + C))', sep = '')
        }
      if(S$Prototroph[i] == T){
            AminoAcidProd <- paste(AminoAcidProd, '+ (A / (A + 2e-03))*', Si, sep = '') #Amino acid production over time
        }
      
    if(BulkSoil == T){
      d_Si <- paste('d_S', i, sep = '') #Begin assigning properties to each species, and building up Glucose and Amino Acid demands
      if(S$LS[i] == 'Copio' & S$Cellulolytic[i] == T &  S$Prototroph[i] == T){
        props <- paste('((1.2*(G/(CKm + G)))*(',Si,'^((Cumax/0.57)+(1-CCP)))) - gamma*',Si, sep = '')
        d_Sidf$d_Si[i] <- paste(props, sep = '')
        GlucDemand <- paste(GlucDemand, '-(1.2*(G/(CKm + G)))*',Si, sep = '')
      }else if(S$LS[i] == 'Copio' & S$Cellulolytic[i] == T &  S$Prototroph[i] == F){
        props <- paste('(((1.2*(G/(CKm + G)))*(1.2*(A/CKm + A)))*(',Si,'^((Cumax/0.57)+(1-CCA)))) - gamma*',Si, sep = '')
        d_Sidf$d_Si[i] <- paste(props, sep = '')
        GlucDemand <- paste(GlucDemand, '-(1.2*(G/(CKm + G)))*',Si, sep = '')
        AminoAcidDemand <- paste(AminoAcidDemand, '-(1.2*(A/CKm + A))*',Si, sep = '')
      }else if(S$LS[i] == 'Copio' & S$Cellulolytic[i] == F &  S$Prototroph[i] == T){
        props <- paste('((1.2*(G/(CKm + G)))*(',Si,'^((Cumax/0.57)+(1-CNP)))) - gamma*',Si, sep = '')
        d_Sidf$d_Si[i] <- paste(props, sep = '')
        GlucDemand <- paste(GlucDemand, '-(1.2*(G/(CKm + G)))*',Si, sep = '')
      }else if(S$LS[i] == 'Copio' & S$Cellulolytic[i] == F &  S$Prototroph[i] == F){
        props <- paste('(((1.2*(G/(CKm + G)))*(1.2*(A/CKm + A)))*(',Si,'^((Cumax/0.57)+(1-CNA)))) - gamma*',Si, sep = '')
        d_Sidf$d_Si[i] <- paste(props, sep = '')
        GlucDemand <- paste(GlucDemand, '-(1.2*(G/(CKm + G)))*',Si, sep = '')
        AminoAcidDemand <- paste(AminoAcidDemand, '-(1.2*(A/CKm + A))*',Si, sep = '')
      }else if(S$LS[i] == 'Oligo' & S$Cellulolytic[i] == T &  S$Prototroph[i] == T){
        props <- paste('((1.2*(G/(OKm + G)))*(',Si,'^((Oumax/0.57)+(1-OCP)))) - gamma*',Si, sep = '')
        d_Sidf$d_Si[i] <- paste(props, sep = '')
        GlucDemand <- paste(GlucDemand, '-(1.2*(G/(OKm + G)))*',Si, sep = '')
      }else if(S$LS[i] == 'Oligo' & S$Cellulolytic[i] == T &  S$Prototroph[i] == F){
        props <- paste('(((1.2*(G/(OKm + G)))*(1.2*(A/OKm + A)))*(',Si,'^((Oumax/0.57)+(1-OCA)))) - gamma*',Si, sep = '')
        d_Sidf$d_Si[i] <- paste(props, sep = '')
        GlucDemand <- paste(GlucDemand, '-(1.2*(G/(OKm + G)))*',Si, sep = '')
        AminoAcidDemand <- paste(AminoAcidDemand, '-(1.2*(A/OKm + A))*',Si, sep = '')
      }else if(S$LS[i] == 'Oligo' & S$Cellulolytic[i] == F &  S$Prototroph[i] == T){
        props <- paste('((1.2*(G/(OKm + G)))*(',Si,'^((Oumax/0.57)+(1-ONP)))) - gamma*',Si, sep = '')
        d_Sidf$d_Si[i] <- paste(props, sep = '')
        GlucDemand <- paste(GlucDemand, '-(1.2*(G/(OKm + G)))*',Si, sep = '')
      }else if(S$LS[i] == 'Oligo' & S$Cellulolytic[i] == F &  S$Prototroph[i] == F){
        props <- paste('(((1.2*(G/(OKm + G)))*(1.2*(A/OKm + A)))*(',Si,'^((Oumax/0.57)+(1-ONA)))) - gamma*',Si, sep = '')
        d_Sidf$d_Si[i] <- paste(props, sep = '')
        GlucDemand <- paste(GlucDemand, '-(1.2*(G/(OKm + G)))*',Si, sep = '')
        AminoAcidDemand <- paste(AminoAcidDemand, '-(1.2*(A/OKm + A))*',Si, sep = '')
      }
    }else{
      Si <- paste('S', i, sep = '')
      d_Si <- paste('d_S', i, sep = '') #Begin assigning properties to each species, and building up Glucose and Amino Acid demands
      if(S$LS[i] == 'Copio' & S$Cellulolytic[i] == T &  S$Prototroph[i] == T){
        props <- paste('((1.2*(G/(CKm + G)))*(',Si,'^((Cumax/0.57)+(1-CCP)))) - gamma*',Si, sep = '')
        d_Sidf$d_Si[i] <- paste(props, sep = '')
        GlucDemand <- paste(GlucDemand, '-(1.2*(G/(CKm + G)))*',Si, sep = '')
      }else if(S$LS[i] == 'Copio' & S$Cellulolytic[i] == T &  S$Prototroph[i] == F){
        props <- paste('(((1.2*(G/(CKm + G))))*(',Si,'^((Cumax/0.57)+(1-CCA)))) - gamma*',Si, sep = '')
        d_Sidf$d_Si[i] <- paste(props, sep = '')
        GlucDemand <- paste(GlucDemand, '-(1.2*(G/(CKm + G)))*',Si, sep = '')
        AminoAcidDemand <- paste(AminoAcidDemand, '-(1.2*(A/CKm + A))*',Si, sep = '')
      }else if(S$LS[i] == 'Copio' & S$Cellulolytic[i] == F &  S$Prototroph[i] == T){
        props <- paste('((1.2*(G/(CKm + G)))*(',Si,'^((Cumax/0.57)+(1-CNP)))) - gamma*',Si, sep = '')
        d_Sidf$d_Si[i] <- paste(props, sep = '')
        GlucDemand <- paste(GlucDemand, '-(1.2*(G/(CKm + G)))*',Si, sep = '')
      }else if(S$LS[i] == 'Copio' & S$Cellulolytic[i] == F &  S$Prototroph[i] == F){
        props <- paste('(((1.2*(G/(CKm + G))))*(',Si,'^((Cumax/0.57)+(1-CNA)))) - gamma*',Si, sep = '')
        d_Sidf$d_Si[i] <- paste(props, sep = '')
        GlucDemand <- paste(GlucDemand, '-(1.2*(G/(CKm + G)))*',Si, sep = '')
        AminoAcidDemand <- paste(AminoAcidDemand, '-(1.2*(A/CKm + A))*',Si, sep = '')
      }else if(S$LS[i] == 'Oligo' & S$Cellulolytic[i] == T &  S$Prototroph[i] == T){
        props <- paste('((1.2*(G/(OKm + G)))*(',Si,'^((Oumax/0.57)+(1-OCP)))) - gamma*',Si, sep = '')
        d_Sidf$d_Si[i] <- paste(props, sep = '')
        GlucDemand <- paste(GlucDemand, '-(1.2*(G/(OKm + G)))*',Si, sep = '')
      }else if(S$LS[i] == 'Oligo' & S$Cellulolytic[i] == T &  S$Prototroph[i] == F){
        props <- paste('(((1.2*(G/(OKm + G))))*(',Si,'^((Oumax/0.57)+(1-OCA)))) - gamma*',Si, sep = '')
        d_Sidf$d_Si[i] <- paste(props, sep = '')
        GlucDemand <- paste(GlucDemand, '-(1.2*(G/(OKm + G)))*',Si, sep = '')
        AminoAcidDemand <- paste(AminoAcidDemand, '-(1.2*(A/OKm + A))*',Si, sep = '')
      }else if(S$LS[i] == 'Oligo' & S$Cellulolytic[i] == F &  S$Prototroph[i] == T){
        props <- paste('((1.2*(G/(OKm + G)))*(',Si,'^((Oumax/0.57)+(1-ONP)))) - gamma*',Si, sep = '')
        d_Sidf$d_Si[i] <- paste(props, sep = '')
        GlucDemand <- paste(GlucDemand, '-(1.2*(G/(OKm + G)))*',Si, sep = '')
      }else if(S$LS[i] == 'Oligo' & S$Cellulolytic[i] == F &  S$Prototroph[i] == F){
        props <- paste('(((1.2*(G/(OKm + G))))*(',Si,'^((Oumax/0.57)+(1-ONA)))) - gamma*',Si, sep = '')
        d_Sidf$d_Si[i] <- paste(props, sep = '')
        GlucDemand <- paste(GlucDemand, '-(1.2*(G/(OKm + G)))*',Si, sep = '')
        AminoAcidDemand <- paste(AminoAcidDemand, '-(1.2*(A/OKm + A))*',Si, sep = '')
      }
    }
    }
  d_C <- paste('C',allCellul, sep = '') #Change in Cellulose over time, controlled by cellulase production
  d_G <- paste(GlucProd, GlucDemand, sep = '') #Compile Glucose production vs use
  d_A <- paste(AminoAcidProd, AminoAcidDemand, sep = '') #Compile Amino Acid production vs use
  return(c(d_C, d_G, d_A, d_Sidf))
}

