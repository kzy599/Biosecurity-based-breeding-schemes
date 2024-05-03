#!~/bin/Rscript
#globe script
for(r in 1:10){

rm(list=setdiff(ls(),"r"))
gc()

#magnitude of GEI
mGEI = 2 # 5 8 represent the genetic correlation 0.2, 0.5, 0.8 respectively

allfiles = list.files()

FounderG0 = paste("GEI",mGEI,"base&g",r,".rda",sep = "")

if(any(allfiles == FounderG0)){load(FounderG0)}else{source("loadfile.R")}

#set parameters
source("parameters.R")
source("package.R")

#output
output <- data.frame(generation=c(0:nGeneration),
                     gv=c(numeric(length = nGeneration+1)),
                     pgv=c(numeric(length = nGeneration+1)),
                     candpgv=c(numeric(length = nGeneration+1)),
                     gv1=c(numeric(length = nGeneration+1)),
                     pgv1=c(numeric(length = nGeneration+1)),
                     candpgv1=c(numeric(length = nGeneration+1)),
                     gxe=c(numeric(length = nGeneration+1)),
                     env1pheno=c(numeric(length = nGeneration+1)),
                     env2pheno=c(numeric(length = nGeneration+1)),
                     env1pgc=c(numeric(length = nGeneration+1)),
                     env2pgc=c(numeric(length = nGeneration+1)),
                     e1e2pc=c(numeric(length = nGeneration+1)),
                     e1pe2g=c(numeric(length = nGeneration+1)),
                     e1ge2p=c(numeric(length = nGeneration+1)),
                     Inbreeding=c(numeric(length = nGeneration+1)),
                     Ne=c(numeric(length = nGeneration+1)),
                     Va=c(numeric(length = nGeneration+1)),
                     Vap1=c(numeric(length = nGeneration+1)),
                     Vap2=c(numeric(length = nGeneration+1)),
                     Inbreeding_plink=c(numeric(length = nGeneration+1)),
                     LD=c(numeric(length = nGeneration+1)),
                     LDscore=c(numeric(length = nGeneration+1)),
                     nCoancestor=c(numeric(length = nGeneration+1)),
                     env1h2=c(numeric(length = nGeneration+1)),
                     env2h2=c(numeric(length = nGeneration+1)),
                     accuracy=c(numeric(length = nGeneration+1)))
#genomic bio-security
if(sch=="bs"){
  source("bs.R")
}

#genomic direct selection
if(sch=="ds"){
  source("ds.R")
}

#traditionl direct selection
if(sch=="tds"){
  source("tds.R")
}

#traditional bio-security selection
if(sch=="tbs"){
  source("tbs.R")
}
}