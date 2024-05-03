pop_select <- selectWithinFam(pop = pop_test,nInd = 10,use = "rand",simParam = SP)
keep<-paste("Animal_",sprintf("%0*d",9,as.numeric(pop_select@id)),sep = "")

pop_LD<- pop_select
pop_Inb <-mergePops(list(pop_founder,pop_select)) 
if(!exists("snpfre_v")){
  snp_012_dt = pullSnpGeno(pop_founder)
  snpfre_v <- apply(snp_012_dt, 2, function(x){
  single_snpfre_s <- sum(x, na.rm = TRUE)/(2*length(x))
  return(single_snpfre_s)
})
}

snp_012_dt = pullSnpGeno(pop_select)

inbreeding_plink = mean(calinbya2(snp_012_dt = snp_012_dt, snpfre = snpfre_v))

ped <-  fread(
  input = "ped.csv",
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)
Pedig <- prePed(ped, keep=keep)
Res   <- pedInbreeding(Pedig)
inbreeding <- mean(Res$Inbr[Res$Indiv %in% keep])


writePlink(pop_LD,baseName = "popLD",simParam = SP)
LDped<- fread(input = "popLD.ped",sep = " ")
LDped[,V3:=0]
LDped[,V4:=0]
fwrite(LDped,file = "popLD.ped",col.names = FALSE,sep = " ")
system("plink --file popLD  --make-bed --chr-set 44 --out popLD")
system("gcta64 --bfile popLD --autosome --ld-score --ld-wind 10000 --ld-rsq-cutoff 0.01 --out popLD")
LDfile<- fread(input="popLD.score.ld",sep = " ")
LD <- mean(LDfile$mean_rsq)
LDscore <- mean(LDfile$ldscore)


if (g == 0) {
  Ne <- 150
} else
  if (g > 0) {
    ped <-  fread(
      input = "ped.csv",
      sep = ",",
      header = TRUE,
      stringsAsFactors = FALSE
    )
    
    pedig <- prePed(ped)
    pKin   <- pedIBD(pedig, keep.only = keep)
    Summary <- summary(pedig)
    id     <- keep
    x      <- Summary[Indiv %in% id]$equiGen
    N      <- length(x)
    n      <- (matrix(x, N, N, byrow = TRUE) + matrix(x, N, N, byrow = FALSE)) / 2
    deltaC <- 1 - (1 - pKin[id, id]) ^ (1 / n)
    Ne   <- 1 / (2 * mean(deltaC))
  }

ped <-  fread(
  input = "ped.csv",
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)
ped_tidy<- visPedigree::tidyped(ped,cand = keep)
output$nCoancestor[g+1]<- length(ped_tidy[Gen==1,Ind])
output$gxe[g+1] <- cor(pop@gv[,1],pop@gv[,2])#varA(pop)[1,2]/(sqrt(varA(pop)[1,1])*sqrt(varA(pop)[2,2]))
output$gv[g+1] <- mean(pop@gv[,2])
output$gv1[g+1] <- mean(pop@gv[,1])
output$Va[g+1] <- varA(pop)[2,2]
output$Ne[g+1] <- Ne
output$Inbreeding[g+1] <- inbreeding
output$Inbreeding_plink[g+1] <- inbreeding_plink
output$LD[g+1] <- LD
output$LDscore[g+1] <- LDscore
if(g==0){
igdt=cali2g2(pop = pop)
vadt=calva(pop=pop)
output$candpgv1[g+1]= mean(pop_founder@gv[,1])
output$candpgv[g+1]= mean(pop_founder@gv[,2])
output$pgv1[g+1]= mean(founderparents@gv[,1])
output$pgv[g+1]= mean(founderparents@gv[,2])
idf=calidf(pop)
famacc=calidf_ebv_gv(pop_test,trait=2)
}else{
igdt=rbind(igdt,cali2g2(pop = pop))
vadt=rbind(vadt,calva(pop=pop))
idf=c(idf,calidf(pop))
famacc=c(famacc,calidf_ebv_gv(pop_test,trait=2))
output$candpgv1[g+1]= mean(candidate@gv[,1])
output$candpgv[g+1]= mean(candidate@gv[,2])
output$pgv1[g+1]= mean(nucleusparents@gv[,1])
output$pgv[g+1]= mean(nucleusparents@gv[,2])
}

output$env1pheno[g+1]=mean(pop@pheno[,1])
output$env2pheno[g+1]=mean(pop@pheno[,2])

output$env1pgc[g+1]=cor(pop@pheno[,1],pop@gv[,1])
output$env2pgc[g+1]=cor(pop@pheno[,2],pop@gv[,2])
output$e1e2pc[g+1]=cor(pop@pheno[,1],pop@pheno[,2])

output$e1pe2g[g+1]=cor(pop@pheno[,1],pop@gv[,2])
output$e1ge2p[g+1]=cor(pop@gv[,1],pop@pheno[,2])

output$env1h2[g+1]=varA(pop)[1,1]/varP(pop)[1,1]
output$env2h2[g+1]=varA(pop)[2,2]/varP(pop)[2,2]

output$vap1[g+1]=varP(pop)[1,1]
output$vap2[g+1]=varP(pop)[2,2]

if(sch%in%c("bs","tbs")){
  output$accuracy[g+1] <- cor(pop_gsb@ebv,pop_gsb@gv[,2])
}else{
  output$accuracy[g+1] <- cor(pop_gsc@ebv,pop_gsc@gv[,2]) 
}
