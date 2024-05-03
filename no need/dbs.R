#direct selection

for (g in 0:nGeneration) {
  
  #population

  pop_test_M <- selectWithinFam(pop = pop,nInd = nProgeny/4,use = "rand",sex = "M",simParam = SP)
  pop_test_F <- selectWithinFam(pop = pop,nInd = nProgeny/4,use = "rand",sex = "F",simParam = SP)
  pop_test <- c(pop_test_F,pop_test_M)

  if(g == 0) rm(.Random.seed) # ensure totally random
  #calculate ebv
  if(exists(x = "test_dt")){
    nucleusparents<- candidate[unique(c(pop@mother,pop@father))]
    test_dt1 <- as.data.table(cbind(paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@id)),sep = ""),
                                    paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@father)),sep = ""),
                                    paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@mother)),sep = ""),
                                    pheno(pop_test)[,1],
                                    pop_test@sex,
                                    paste("G",g,sep = ""),
                                    pop_test@gv[,1]))
    colnames(test_dt1) <- c("AnimalID","SireID","DamID","pheno","Sex","Generation","gv")
    FamilyID<- paste("Family_",substr(test_dt1$DamID,start = 8,stop = 16),"_",substr(test_dt1$SireID,start = 8,stop =16 ),sep = "")
    test_dt1[,Family:=FamilyID]
    if(!is.null(test_dt$ebv)){
      test_dt[,ebv:=NULL]
    }
    test_dt<- rbind(test_dt,test_dt1)
  }else{
    founderparents<- pop_founder[unique(c(pop@mother,pop@father))]
    test_dt1 <- as.data.table(cbind(paste("Animal_",sprintf("%0*d",9,as.numeric(founderparents@id)),sep = ""),
                                    NA,
                                    NA,
                                    NA,
                                    founderparents@sex,
                                    "founder",
                                    founderparents@gv[,1]))
    colnames(test_dt1) <- c("AnimalID","SireID","DamID","pheno","Sex","Generation","gv")
    test_dt1[,Family:=NA]
    test_dt <- as.data.table(cbind(paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@id)),sep = ""),
                                   paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@father)),sep = ""),
                                   paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@mother)),sep = ""),
                                   pheno(pop_test)[,1],
                                   pop_test@sex,
                                    paste("G",g,sep = ""),
                                   pop_test@gv[,1]))
    colnames(test_dt) <- c("AnimalID","SireID","DamID","pheno","Sex","Generation","gv")
    FamilyID<- paste("Family_",substr(test_dt$DamID,start = 8,stop = 16),"_",substr(test_dt$SireID,start = 8,stop =16 ),sep = "")
    test_dt[,Family:=FamilyID]
    test_dt<- rbind(test_dt1,test_dt)
  }
  
  dt <- test_dt
  setorder(dt,AnimalID)
  ped <- dt[,1:3]
  ped<- visPedigree::tidyped(ped = ped,cand = dt$AnimalID)
  
  fwrite(dt,file = "pheno.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "0")
  fwrite(ped[,1:3],file = "ped.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "0")
  
  if(!exists("p_geno")){
  p_geno <- pullSnpGeno(pop_founder[unique(c(pop@mother,pop@father))])
  rownames(p_geno)<- paste("Animal_",sprintf("%0*d",9,as.numeric(rownames(p_geno))),sep = "") 
}else{
  nucleusparentsgeno <- pullSnpGeno(nucleusparents)
  rownames(nucleusparentsgeno)<- paste("Animal_",sprintf("%0*d",9,as.numeric(rownames(nucleusparentsgeno))),sep = "") 
  p_geno <- rbind(p_geno,nucleusparentsgeno)
}
  
  pheno_file_s <- "pheno.csv"
  ped_file_s <- "ped.csv"
 
  pheno_dt <-
    fread(
      pheno_file_s,
      sep = ",",
      header = TRUE,
      stringsAsFactors = FALSE,
      na.strings = "NA"
    )
  ped_dt <-
    fread(
      ped_file_s,
      sep = ",",
      header = TRUE,
      stringsAsFactors = FALSE,
      na.strings = "NA"
    )
  #################
  fwrite(
    pheno_dt,
    file = "pheno.txt",
    na = "0",
    append = FALSE,
    sep = " ",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  fwrite(
    ped_dt[,1:3],
    file = "ped.txt",
    append = FALSE,
    sep = " ",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE,
    na="0"
  )
  
    source("blup_ped.R")
  
  
  pop_test@ebv<- testfunction(popname = pop_test,test1 = dt)
  
  
  candfam_reference<- selectFam(pop = pop_test,nFam = ncf,use = "ebv",simParam = SP)
  
  
  #genotype
  pop_gsc_male <- selectWithinFam(pop = candfam_reference,nInd = round(50*(1/3)),use = "ebv",simParam = SP,sex="M")
  pop_gsc_female <- selectWithinFam(pop = candfam_reference,nInd = round(50*(2/3)),use = "ebv",simParam = SP,sex="F")
  pop_gsc = c(pop_gsc_female,pop_gsc_male)
  
  conferencegeno<- pullSnpGeno(pop_gsc,snpChip = 1,simParam = SP)
  rownames(conferencegeno) <- paste("Animal_",sprintf("%0*d",9,as.numeric(rownames(conferencegeno))),sep = "")
  
  if(!exists("geno")){
  geno = conferencegeno
}else{
 geno <- rbind(geno,conferencegeno)
}
  
  if(g>3){
    genoname<- dt[Generation%in%paste("G",c((g-3):g),sep = ""),AnimalID]
    geno<- geno[rownames(geno)%in%genoname,]
  }
  pgenoneed= p_geno[!rownames(p_geno)%in%rownames(geno),]
  genodt<- rbind(pgenoneed,geno)
  genodt<- as.data.table(genodt,keep.rownames = TRUE)
  rn<- genodt$rn
  genodt<- genodt[,-1]
  genodt <- apply(genodt,1,setrowcollapse)
  genodt <-data.table(rn,genodt)
  fwrite(genodt,file = "geno_selectedparents.txt",sep = " ",col.names = FALSE, row.names = FALSE,quote = FALSE)
  
  
  source("blup.R")
  
  
  pop_gsc@ebv<- testfunction(popname = pop_gsc,test1 = dt)
  
  source("population_parameters_calculation_dbs.R")
  
  candidate_female<-  selectWithinFam(pop_gsc,nInd = 4,use = "ebv",sex="F",simParam = SP)
  candidate_male<-  selectWithinFam(pop_gsc,nInd = 2,use = "ebv",sex="M",simParam = SP)
  candidate=c(candidate_female,candidate_male)
  if(g<nGeneration){
    pop <- ocs(
      pop = candidate,
      nCrosses = nfamily,
      nProgenyPerCross = nProgeny,
      nFemalesMax = 100,
      nMalesMax = 50,
      equalizeFemaleContributions =TRUE,
      equalizeMaleContributions = TRUE,
      targetDegree = 45,
      use = "ebv"
    )  
  }
  
  
}
famdt=data.table(idf,famacc)
output=cbind(output,famdt,igdt,vadt)
fwrite(output,file = paste(bsname,r,".csv",sep = ""),sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE)
