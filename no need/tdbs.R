#tds
for (g in 0:nGeneration) {
  
  #population
  pop_test_M <- selectWithinFam(pop = pop,nInd = nProgeny/4,use = "rand",sex = "M",simParam = SP)
  pop_test_F <- selectWithinFam(pop = pop,nInd = nProgeny/4,use = "rand",sex = "F",simParam = SP)
  pop_test <- mergePops(list(pop_test_F,pop_test_M))

  if(g == 0) rm(.Random.seed) # ensure totally random
  
  #calculate ebv
  if(exists(x = "test_dt")){
    nucleusparents<- candidate[unique(c(pop@mother,pop@father))]
    test_dt1 <- as.data.table(cbind(paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@id)),sep = ""),
                                    paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@father)),sep = ""),
                                    paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@mother)),sep = ""),
                                    pheno(pop_test)[,1],
                                    pop_test@sex,
                                    pop_test@gv[,1]))
    colnames(test_dt1) <- c("AnimalID","SireID","DamID","pheno","Sex","gv")
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
                                    founderparents@gv[,1]))
    colnames(test_dt1) <- c("AnimalID","SireID","DamID","pheno","Sex","gv")
    test_dt1[,Family:=NA]
    test_dt <- as.data.table(cbind(paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@id)),sep = ""),
                                   paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@father)),sep = ""),
                                   paste("Animal_",sprintf("%0*d",9,as.numeric(pop_test@mother)),sep = ""),
                                   pheno(pop_test)[,1],
                                   pop_test@sex,
                                   pop_test@gv[,1]))
    colnames(test_dt) <- c("AnimalID","SireID","DamID","pheno","Sex","gv")
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
    7.75,
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
    5.39
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
  
  pop_test@ebv<- testfunction(popname = pop_test,test1 = dt)
  
  
  candfam_reference<- selectFam(pop = pop_test,nFam = ncf,use = "ebv",simParam = SP)
  
  
  #genotype
  if(gsc=="TOP"){
    pop_gsc_male <- selectWithinFam(pop = candfam_reference,nInd = round(ngi*(1/3)),use = "ebv",simParam = SP,sex="M")
    pop_gsc_female <- selectWithinFam(pop = candfam_reference,nInd = round(ngi*(2/3)),use = "ebv",simParam = SP,sex="F")
    pop_gsc = c(pop_gsc_female,pop_gsc_male)
  }else if(gsc=="TB"){
    pop_gsc1 <- selectWithinFam(pop = candfam_reference,nInd = ngi/2,trait = 1,use = "ebv",simParam = SP)
    pop_gsc2 <- selectWithinFam(pop = candfam_reference,nInd = ngi/2,trait = 1,use = "ebv",selectTop = FALSE,simParam = SP)
    pop_gsc<- mergePops(list(pop_gsc1,pop_gsc2))
  }else if(gsc=="RAND"){
    pop_gsc <- selectWithinFam(pop = candfam_reference,nInd = ngi,trait = 1,use = "ebv",simParam = SP)
  }
  
  source("population_parameters_calculation_dbs.R")
  
  candidate_female<-  selectWithinFam(pop_gsc,nInd = 4,use = "ebv",sex="F",simParam = SP)
  candidate_male<-  selectWithinFam(pop_gsc,nInd = 2,use = "ebv",sex="M",simParam = SP)
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