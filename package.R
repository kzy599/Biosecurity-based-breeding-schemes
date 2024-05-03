library(AlphaSimR)
library(data.table)
library(visPedigree)
library(optiSel)
library(nadiv)
library(sampling)
library(AlphaMME)
library(AlphaLearn)
setrowcollapse <- function(x){
  return(paste(x,sep="",collapse = ""))
}
calinbya2 <- function(snp_012_dt,snpfre){
  calformula <- function(x, snpfre){
    y <- mean((x^2 - ((1 + 2 * snpfre) * x) + 2 * (snpfre^2)) / (2 * snpfre * (1 - snpfre)), na.rm = TRUE)
    return(y)
  }
  inb <- apply(snp_012_dt, 1, calformula,snpfre = snpfre)
  return(inb)
}

abstractebv<- function(popname,dtebv){
  id_temp <- paste("Animal_",sprintf("%0*d",9,as.numeric(popname@id)),sep = "")
  popebv = dtebv$ebv[match(id_temp,dtebv$AnimalID)]
  popebv<- matrix(as.numeric(popebv),nrow = popname@nInd,ncol = 1)
  return(popebv)
}
  cali2g2<- function(pop){
    va1<- varA(pop)[1,1]
    va2<- varA(pop)[2,2]
    cov12 <- varA(pop)[1,2]
    g2<- ((va1+va2)/2+cov12)/2
    i2 <- g2-cov12
    x<- data.table(g2,i2)
    return(x)
  }
  calva <- function(pop){
    va1<- varA(pop)[1,1]
    va2<- varA(pop)[2,2]
    cov12 <- varA(pop)[1,2]
    x<- data.table(va1,va2,cov12)
    return(x)
  }

#calculate the accuracy of family rank prediction
calidf_ebv_gv <- function(pop,trait){
popfam1<- selectFam(pop,nFam = 50,trait = trait,use = "gv")
fam1<- unique(paste(popfam1@mother,popfam1@father))

popfam2<- selectFam(pop,nFam = 50,trait = 1,use = "ebv",simParam = SP)
fam2<- unique(paste(popfam2@mother,popfam2@father))
x<-sum(fam1 %in% fam2)
return(x)
}

#calculate the number of identical top 50 family between env1 and env2
calidf <- function(pop){
  popfam1<- selectFam(pop,nFam = 50,trait = 2,use = "gv")
  fam1<- unique(paste(popfam1@mother,popfam1@father))
  
  popfam2<- selectFam(pop,nFam = 50,trait = 1,use = "gv",simParam = SP)
  fam2<- unique(paste(popfam2@mother,popfam2@father))
  x<- sum(fam1 %in% fam2)
  return(x)
}

#abstractebv = function(pop,hiebv){
 # colnames(hiebv)[2] = "hib.HA"
  #ebv = hiebv$hib.HA[match(pop@id,hiebv$ID)]
  #ebv = matrix(ebv,ncol = 1)
  #return(ebv)
#}

#abstractebv_blupf90 = function(pop,dt){
 # reqid = paste("Animal_",sprintf("%0*d",9,as.numeric(pop@id)),sep = "")
 # popebv = dt$ebv[match(reqid,dt$AnimalID)]
 # popebv = matrix(popebv,ncol = 1)
 # return(popebv)
#}


#abstract additive effect of individuals from the result of hiblup for pblup
testfunction_hiblup_ped<- function(popname,test1){
    test2<- as.integer(popname@id)
    test2<- as.data.table(test2)
    colnames(test2) <- "ID"
    test1<- test1[ID%in%test2$ID,]
    test2<- merge.data.table(x=test2,y=test1,by = "ID",sort = FALSE)
    test2 <- matrix(test2$hib.PA,nrow = popname@nInd,ncol = 1)
    return(test2)
  }

  #Calculate mean fre and square
  #The finalmat have 1001 lines that first line contain the minor frequence of each snp 
  calfreandmat=function(pop){
    mat = pullSnpGeno(pop)
    fre=apply(mat,2,function(x){
      y=sum(x)/(2*length(x))
      return(y)
    })
    fre[fre>0.5]=1-fre[fre>0.5]
    fremat=matrix(rep(fre,nrow(mat)),ncol=length(fre),nrow=nrow(mat),byrow=TRUE)
    zmat=mat-fremat
    finalmat=rbind(fre,zmat)
    return(finalmat)
  }
  calsnploci = function(finalmat){
    fre=finalmat[1,]
    lociname = colnames(finalmat)[fre>0.05]
    return(lociname)
  }

  calfre=function(finalmat){
    y=mean(finalmat[1,])
    return(y)
  }

  calsquare=function(finalmat){
    z=finalmat[-1,]
    zz= z %*% t(z)
    sq=mean(diag(zz))
    return(sq)
  }
  calhet=function(pop){
    z=pullSnpGeno(pop)
    het=apply(z,1,function(x){
      y=sum(x*(2-x))/length(x)
      return(y)
    })
    meanhet=mean(het)
    return(meanhet)
  }

#Funtion about creat file of plink is not used in this simulation
  #make .ped file required by plink
  makeped = function(pop1,pop2,calparents,pop_parents){
    z1=pullSnpGeno(pop1)
    fre1=apply(z1,2,function(x){
      y=sum(x)/(2*length(x))
      return(y)
    })
    fre1[fre1>0.5]=1-fre1[fre1>0.5]
    snploci1 = names(fre1)[fre1>0.05]
    
    if(calparents==TRUE){
      pop2 = mergePops(list(pop_parents,pop2))
    }
    z2=pullSnpGeno(pop2)
    fre2=apply(z2,2,function(x){
      y=sum(x)/(2*length(x))
      return(y)
    })
    fre2[fre2>0.5]=1-fre2[fre2>0.5]
    snploci2 = names(fre2)[fre2>0.05]

    snplocifinal = unique(c(snploci1,snploci2))

    if(calparents==TRUE){
     pop = mergePops(list(pop1,pop2))
    }else{
     pop = mergePops(list(pop_parents,pop1,pop2))
    }
    z=pullSnpGeno(pop)
    z=z[,colnames(z)%in%snplocifinal]
    
    fre=apply(z,2,function(x){
      y=sum(x)/(2*length(x))
      return(y)
    })
    
    frename=names(fre)[fre>0.5]
    #idped=rownames(z)
    
    z[,frename] = (2-z[,frename]) 
    
    z[z==2]=22
    z[z==1]=12
    z[z==0]=11

    z = as.data.table(cbind(pop@id,z))
    return(z)
  }

  #make .map file required by plink 
  makemap = function(...){
  mapck = getSnpMap()
  mapck = mapck[colnames(genodt)[-1],]
  mapck$id = mapck$chr
  mapck$chr = rownames(mapck)
  mapck$site = mapck$pos
  mapck$pos = rep(0,nrow(mapck))
  return(mapck)
  }

source("ocs.R")
source("ocsped.R")