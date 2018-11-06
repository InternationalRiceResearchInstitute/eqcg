#####################
# getEqcgs
#' Function to calculate the number eqivalent complete generations from a pedigree
#'
#' Parameters
#' @param ped is an pedigree dataframe consiting of the columns label, sire, and dam
#' @param targ is a vector of ids for the individuals of interest
#' 
#' @return a vector of eqcg values for the individuals of interest
#####################

getEqcgs<- function(ped, targ){

  #get parents function
  getPrnts<- function(gid, ped){
    ix<- match(gid, ped@label)
    sr0<- ped@sire[ix]
    dm0<- ped@dam[ix]
    sr<- ifelse(is.na(sr0),0,sr0)
    dm<- ifelse(is.na(dm0),0,  dm0)
    while(sr==dm & (sr+dm) >0){
      ix<- sr
      sr0<- ped@sire[ix]
      dm0<- ped@dam[ix]
      sr<- ifelse(is.na(sr0),0, sr0)
      dm<- ifelse(is.na(dm0),0, dm0)
    }
    if(dm>0){
      dm<- ped@label[dm]
    }
    if(sr>0){
      sr<- ped@label[sr]
    }
    return(as.numeric(c(dm, sr)))
  }
  
  
  #Eqcg
  calcEqcg<- function(GID, pd){
    if(GID==0){
      start<- 0
    }else{
      cyc<-0
      p<- getPrnts(GID, ped=pd)
      p<- p[which(p!=0)]
      if(length(p)>0){
        cyc<- cyc+1
        inc<- sum(0.5^cyc * matrix(1, length(p)))
        start<- inc
        while(inc>0){
          p<- unlist(sapply(p, getPrnts, pd))
          cyc<- cyc+1
          p<- p[which(p!=0)]
          inc<- 0
          if(length(p)>0){
            inc<- sum( 0.5^cyc * matrix(1, length(p)) )
          }
          start<- start+inc
        }
      }else{
        start<- 0
      }
    }
    return(start)
  }
  
  
  #fast Eqcg
  fastEqcg<- function(eqcgs){
    eqcg<- mean(c(eqcgs[1], eqcgs[2]))+1
    return(eqcg)
  }
  
  #Prepare pedigree file
  require(pedigreemm)
  require(nadiv)
  colnames(ped)<- c('label', 'sire', 'dam')
  ped<- nadiv::prepPed(ped)
  ped<- nadiv::prunePed(ped, phenotyped=targ)
  ped<- pedigreemm::editPed(ped$sire, ped$dam, ped$label)
  P3<-  pedigreemm::pedigree(ped$sire, ped$dam,ped$label)
  
  #Loop to use fast EqCg in an ordered pedigree
  ped<- data.frame(ped, eqcg=NA)
  for(i in 1:nrow(ped)){
    pars<- getPrnts(ped[i,1], P3)
    if(0 %in% pars){
      ped[i,'eqcg']<- calcEqcg(ped[i,1], P3)
    }else{
      eqcgPar<- ped[match(pars, ped[,1]), 'eqcg']
      if(NA %in% eqcgPar){
        ped[i,'eqcg']<- calcEqcg(ped[i,1], P3)
      }else{
        ped[i,'eqcg']<- fastEqcg(eqcgPar)
      }
    }
  }  
return(ped[match(targ, ped[,1]),'eqcg'])
}


