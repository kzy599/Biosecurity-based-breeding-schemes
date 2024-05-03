#traditional bio-security

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
                                    nucleusparents@gv[,2]))
    colnames(test_dtp) <- c("AnimalID","SireID","DamID","pheno","Sex","gv")
    FamilyID<- paste("Family_",substr(test_dtp$DamID,start = 8,stop = 16),"_",substr(test_dtp$SireID,start = 8,stop =16 ),sep = "")
    test_dtp[,Family:=FamilyID]
    test_dt1 <- as.data.table(cbind(paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@id)),sep = ""),
                                    paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@father)),sep = ""),
                                    paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@mother)),sep = ""),
                                    pheno(pop_test)[,2],
                                    pop_test@sex,
                                    pop_test@gv[,2]))
    colnames(test_dt1) <- c("AnimalID","SireID","DamID","pheno","Sex","gv")
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
                                    founderparents@gv[,2]))
    colnames(test_dt1) <- c("AnimalID","SireID","DamID","pheno","Sex","gv")
    test_dt1[,Family:=NA]
    test_dt <- as.data.table(cbind(paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@id)),sep = ""),
                                   paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@father)),sep = ""),
                                   paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@mother)),sep = ""),
                                   pheno(pop_test)[,2],
                                   pop_test@sex,
                                   pop_test@gv[,2]))
    colnames(test_dt) <- c("AnimalID","SireID","DamID","pheno","Sex","gv")
    FamilyID<- paste("Family_",substr(test_dt$DamID,start = 8,stop = 16),"_",substr(test_dt$SireID,start = 8,stop =16 ),sep = "")
    test_dt[,Family:=FamilyID]
    test_dt<- rbind(test_dt1,test_dt)
  }
  
  breed_dt <- as.data.table(cbind(paste("Animal_",sprintf("%0*d",9,as.numeric(pop_breed@id)),sep = ""),
                                  paste("Animal_",sprintf("%0*d",9,as.numeric(pop_breed@father)),sep = ""),
                                  paste("Animal_",sprintf("%0*d",9,as.numeric(pop_breed@mother)),sep = ""),
                                  NA,
                                  pop_breed@sex,
                                  pop_breed@gv[,2]))
  colnames(breed_dt) <- c("AnimalID","SireID","DamID","pheno","Sex","gv")
  FamilyID<- paste("Family_",substr(breed_dt$DamID,start = 8,stop = 16),"_",substr(breed_dt$SireID,start = 8,stop =16 ),sep = "")
  breed_dt[,Family:=FamilyID]
  
  dt <- rbind(test_dt,breed_dt)
  setorder(dt,AnimalID)
  ped <- dt[,1:3]
  ped<- visPedigree::tidyped(ped = ped,cand = dt$AnimalID)
  
  fwrite(dt,file = "pheno.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "0")
  fwrite(ped[,1:3],file = "ped.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "0")
  
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
  parameter_txt_v <- c(
    # 指定数据文件名称，默认名称，不能修改
    "DATAFILE",
    "pheno.txt",
    # *待分析性状所在列位置
    "TRAITS",
    4,
    # 追加到输出文件中的列位置，可以为空
    "FIELDS_PASSED TO OUTPUT",
    "",
    # 性状加权，可以为空  
    "WEIGHT(S)",
    "",
    # 残差方差
    "RESIDUAL_VARIANCE",
    varP(pop)[2,2] - varA(pop)[2,2],
    "EFFECT",
    # *个体编号列位置，crossclassified, 字符类型
    "1 cross alpha",
    # 随机效应
    "RANDOM",
    # 个体动物模型
    "animal",
    # 定义系谱文件名称，默认名称，不能修改
    "FILE",
    "ped.txt",
    # *定义个体、父本和母本在系谱文件中的列位置
    "FILE_POS",
    "1 2 3 0 0",
    # 系谱深度，0表示追溯到奠基者世代
    "PED_DEPTH",
    0,
    # 加性遗传方差
    "(CO)VARIANCES",
    varA(pop)[2,2]
  )
  # 生成参数文件renumf90.par，utf-8格式
  parameter_con <- file("renumf90.par")
  writeLines(text = paste(parameter_txt_v, collapse = "\n"),
             con = parameter_con)
  close(parameter_con)
  
  # 调用renumf90程序
  system2("renumf90", args = "renumf90.par", stdout="renumf90.log", wait = TRUE)
  
  system2("blupf90", args = "renf90.par", stdout="blupf90.log",  wait = TRUE)
  solution<- fread(input="solutions",sep = ' ')
  setorder(solution,level)
  rename<- fread(input = "renadd01.ped",sep = " ")
  trueoder<- rename[,.(V1,V10)]
  setorder(trueoder,V1)
  ebv<- solution[`trait/effect`==1,solution]
  trueorder<- cbind(trueoder,ebv)
  setorder(trueorder,V10)
  trueorder<- as.data.table(trueorder)
  dt[,ebv:=trueorder[,ebv]]
  
  pop_test@ebv<- abstractebv(popname = pop_test,dtebv = dt)

  pop_breed@ebv = abstractebv(popname = pop_breed,dtebv = dt)
  
  candfam_reference<- selectFam(pop = pop_test,nFam = ncf,use = "ebv",simParam = SP)
  
  #candfam_breed<- selectFam(pop = pop_breed,nFam = ncf,use = "ebv",simParam = SP)
  
  #fam_cand= unique(paste(candfam_breed@mother,candfam_breed@father,sep = ""))
  #fam_refe= unique(paste(candfam_reference@mother,candfam_reference@father,sep = ""))
  
  candfam_selection<- pop_breed[pop_breed@mother%in%unique(candfam_reference@mother)]
  
  
  
  if(gsb=="TOP"){
    pop_gsb_male <- selectWithinFam(pop = candfam_selection,nInd = round(ngi*(1/3)),trait = 1,use = "pheno",simParam = SP,sex="M")
    pop_gsb_female <- selectWithinFam(pop = candfam_selection,nInd = round(ngi*(2/3)),trait = 1,use = "pheno",simParam = SP,sex="F")
    pop_gsb = c(pop_gsb_female,pop_gsb_male)
  }else if(gsb=="TB"){
    pop_gsb1 <- selectWithinFam(pop = candfam_selection,nInd = ngi/2,trait = 1,use = "pheno",simParam = SP)
    pop_gsb2 <- selectWithinFam(pop = candfam_selection,nInd = ngi/2,trait = 1,use = "pheno",selectTop = FALSE,simParam = SP)
    pop_gsb<- mergePops(list(pop_gsb1,pop_gsb2))
  }else if(gsb=="RAND"){
    pop_gsb <- selectWithinFam(pop = candfam_selection,nInd = ngi,trait = 1,use = "rand",simParam = SP)
  }

  
  #using phenotype to account within family variance
  #using ebv to account between family variance
  stdphe = function(candfam_selection,pop){

  #phenotype in NE, EBV in CE  
  candt = data.table(id = candfam_selection@id,phe = pheno(candfam_selection)[,1],ebv = candfam_selection@ebv,fam = paste(candfam_selection@mother,candfam_selection@father,sep = ""))
  
  colnames(candt) = c("id","phe","ebv","fam")
  
  #Standard normalization of both value for all candidate families
  #there was no environemnt effect for different family
  candt$phe = (candt$phe - mean(candt$phe))/sd(candt$phe)
    
  candt$ebv = (candt$ebv - mean(candt$ebv))/sd(candt$ebv)

  famns = unique(candt$fam)
    
  #minus the family mean from the phenotype, using deviation to account within family variance
  for(std in 1:length(famns)){
   nedphe = candt[fam==famns[std],phe]
   nedphe = nedphe - mean(nedphe)
   candt[fam==famns[std],phe:=nedphe]
  } 

  #plus the family EBV to account among family variance
  #match the value for the selection candidates
  stdebv = candt$ebv[match(pop@id,candt$id)]
  stdphe = candt$phe[match(pop@id,candt$id)]
  result = matrix((stdebv+stdphe),ncol = 1)
  
  #return the value
  return(result)
  }
  
  pop_gsb@ebv = stdphe(candfam_selection = candfam_selection, pop = pop_gsb)
  
  source("population_parameters_calculation.R")
  
  candidate_female<-  selectWithinFam(pop_gsb,nInd = 4,use = "ebv",sex="F",simParam = SP)
  candidate_male<-  selectWithinFam(pop_gsb,nInd = 2,use = "ebv",sex="M",simParam = SP)
  candidate=c(candidate_female,candidate_male)
  if(g<nGeneration){
    pop <- ocsBLUP(
      pop = candidate,
      pedmore = ped,
      pedkf = FALSE,
      number = 9,
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
