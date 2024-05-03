#bio-security

for (g in 0:nGeneration) {
  
#population
pop_test_M <- selectWithinFam(pop = pop,nInd = nProgeny/4,use = "rand",sex = "M",simParam = SP)
pop_test_F <- selectWithinFam(pop = pop,nInd = nProgeny/4,use = "rand",sex = "F",simParam = SP)
pop_test <- c(pop_test_F,pop_test_M)
pop_breed <- pop[!pop@id%in%pop_test@id]

if(g == 0) rm(.Random.seed) # ensure totally random

#calculate ebv
if(exists(x = "test_dt")){
  nucleusparents<- candidate[unique(c(pop@mother,pop@father))]
  test_dtp <- as.data.table(cbind(paste("Animal_",sprintf("%0*d",9,as.numeric(nucleusparents@id)),sep = ""),
                                  paste("Animal_",sprintf("%0*d",9,as.numeric(nucleusparents@father)),sep = ""),
                                  paste("Animal_",sprintf("%0*d",9,as.numeric(nucleusparents@mother)),sep = ""),
                                  NA,
                                  nucleusparents@sex,
                                  paste("G",(g-1),sep = ""),
                                  nucleusparents@gv[,2]))
  colnames(test_dtp) <- c("AnimalID","SireID","DamID","pheno","Sex","Generation","gv")
  FamilyID<- paste("Family_",substr(test_dtp$DamID,start = 8,stop = 16),"_",substr(test_dtp$SireID,start = 8,stop =16 ),sep = "")
  test_dtp[,Family:=FamilyID]
 test_dt1 <- as.data.table(cbind(paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@id)),sep = ""),
                                 paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@father)),sep = ""),
                                 paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@mother)),sep = ""),
                                 pheno(pop_test)[,2],
                                 pop_test@sex,
                                 paste("G",g,sep = ""),
                                 pop_test@gv[,2]))
 colnames(test_dt1) <- c("AnimalID","SireID","DamID","pheno","Sex","Generation","gv")
 FamilyID<- paste("Family_",substr(test_dt1$DamID,start = 8,stop = 16),"_",substr(test_dt1$SireID,start = 8,stop =16 ),sep = "")
 test_dt1[,Family:=FamilyID]
 test_dt<- rbind(test_dt,test_dtp,test_dt1)
}else{
  founderparents<- pop_founder[unique(c(pop@mother,pop@father))]
  test_dt1 <- as.data.table(cbind(paste("Animal_",sprintf("%0*d",9,as.numeric(founderparents@id)),sep = ""),
                                 NA,
                                 NA,
                                 NA,
                                 founderparents@sex,
                                 "founder",
                                 founderparents@gv[,2]))
  colnames(test_dt1) <- c("AnimalID","SireID","DamID","pheno","Sex","Generation","gv")
  test_dt1[,Family:=NA]
  test_dt <- as.data.table(cbind(paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@id)),sep = ""),
                                 paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@father)),sep = ""),
                                 paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@mother)),sep = ""),
                                 pheno(pop_test)[,2],
                                 pop_test@sex,
                                 paste("G",g,sep = ""),
                                 pop_test@gv[,2]))
  colnames(test_dt) <- c("AnimalID","SireID","DamID","pheno","Sex","Generation","gv")
  FamilyID<- paste("Family_",substr(test_dt$DamID,start = 8,stop = 16),"_",substr(test_dt$SireID,start = 8,stop =16 ),sep = "")
  test_dt[,Family:=FamilyID]
  test_dt<- rbind(test_dt1,test_dt)
}

breed_dt <- as.data.table(cbind(paste("Animal_",sprintf("%0*d",9,as.numeric(pop_breed@id)),sep = ""),
                            paste("Animal_",sprintf("%0*d",9,as.numeric(pop_breed@father)),sep = ""),
                            paste("Animal_",sprintf("%0*d",9,as.numeric(pop_breed@mother)),sep = ""),
                            NA,
                            pop_breed@sex,
                            paste("G",g,sep = ""),
                            pop_breed@gv[,2]))
colnames(breed_dt) <- c("AnimalID","SireID","DamID","pheno","Sex","Generation","gv")
FamilyID<- paste("Family_",substr(breed_dt$DamID,start = 8,stop = 16),"_",substr(breed_dt$SireID,start = 8,stop =16 ),sep = "")
breed_dt[,Family:=FamilyID]

dt <- rbind(test_dt,breed_dt)
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
snp_file_s <- "geno_selectedparents.txt"
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

pop_test@ebv<- abstractebv(popname = pop_test,dtebv = dt)


candfam_reference<- selectFam(pop = pop_test,nFam = ncf,use = "ebv",simParam = SP)

candfam_selection<- pop_breed[pop_breed@mother%in%unique(candfam_reference@mother)]


#genotype
#genotype
if(gsc=="TOP"){
  pop_gsc <- selectWithinFam(pop = candfam_reference,nInd = 50,use = "ebv",simParam = SP)
}else if(gsc=="TB"){
  pop_gsc1 <- selectWithinFam(pop = candfam_reference,nInd = 25,use = "ebv",simParam = SP)
  pop_gsc2 <- selectWithinFam(pop = candfam_reference[!candfam_reference@id%in%pop_gsc1@id],nInd = 25,use = "ebv",selectTop = FALSE,simParam = SP)
  pop_gsc<- c(pop_gsc1,pop_gsc2)
}else if(gsc=="RAND"){
  pop_gsc <- selectWithinFam(pop = candfam_reference,nInd = 50,trait = 2,use = "rand",simParam = SP)
}else if(gsc=="TR"){
  pop_gsc1 = selectWithinFam(pop = candfam_reference,nInd = 25,use = "ebv",simParam = SP)
  pop_gsc2 = selectWithinFam(pop = candfam_reference[!candfam_reference@id%in%pop_gsc1@id],nInd = 25,trait = 2,use = "rand",simParam = SP)
  pop_gsc<- c(pop_gsc1,pop_gsc2)
}

pop_gsb_male <- selectWithinFam(pop = candfam_selection,nInd = round(ngi*(1/3)),trait = 1,use = "pheno",simParam = SP,sex="M")

pop_gsb_female <- selectWithinFam(pop = candfam_selection,nInd = round(ngi*(1/3)),trait = 1,use = "pheno",simParam = SP,sex="F")

pop_gsb = c(pop_gsb_female,pop_gsb_male)

conferencegeno<- pullSnpGeno(pop_gsc,snpChip = 1,simParam = SP)
rownames(conferencegeno) <- paste("Animal_",sprintf("%0*d",9,as.numeric(rownames(conferencegeno))),sep = "")

if(!exists("geno")){
  geno = conferencegeno
}else{
 geno <- rbind(geno,conferencegeno)
}

if(g>3){
  genoname<- dt[Generation%in%paste("G",c((g-3):(g)),sep = ""),AnimalID]
  geno<- geno[rownames(geno)%in%genoname,]
}

selectiongeno <- pullSnpGeno(pop_gsb,snpChip = 1,simParam = SP)
rownames(selectiongeno) <- paste("Animal_",sprintf("%0*d",9,as.numeric(rownames(selectiongeno))),sep = "")

genodt<- rbind(p_geno,geno,selectiongeno)

genodt = genodt[rownames(genodt)%in%dt$AnimalID,]

genodt<- as.data.table(genodt,keep.rownames = TRUE)

rn<- genodt$rn
genodt<- genodt[,-1]
genodt <- apply(genodt,1,setrowcollapse)
genodt <-data.table(rn,genodt)
fwrite(genodt,file = "geno_selectedparents.txt",sep = " ",col.names = FALSE, row.names = FALSE,quote = FALSE)

source("blup.R")

pop_gsb@ebv<- abstractebv(popname = pop_gsb,dtebv = dt)

source("population_parameters_calculation.R")

candidate_female<- selectWithinFam(pop_gsb,nInd = 4,use = "ebv",sex="F",simParam = SP)
candidate_male<-  selectWithinFam(pop_gsb,nInd = 2,use = "ebv",sex="M",simParam = SP)
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
