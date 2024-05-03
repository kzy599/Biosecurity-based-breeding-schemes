parameter_txt_v <- c(
  # the name of file of pheno
  "DATAFILE",
  "pheno.txt",
  #the position of trait in the file of pheno
  "TRAITS",
  4,
  #
  "FIELDS_PASSED TO OUTPUT",
  "",
  #
  "WEIGHT(S)",
  "",
  #
  "RESIDUAL_VARIANCE",
  varP(pop)[2,2] - varA(pop)[2,2],
  "EFFECT",
  # random effect, id
  "1 cross alpha",
  #
  "RANDOM",
  # animal model
  "animal",
  # name of the file of ped
  "FILE",
  "ped.txt",
  # the position of id sir dam
  "FILE_POS",
  "1 2 3 0 0",
  "SNP_FILE",
  "geno_selectedparents.txt",
  # ped depth
  "PED_DEPTH",
  0,
  # additive genetic variance
  "(CO)VARIANCES",
  varA(pop)[2,2],
  "OPTION minfreq 0.05",
  "OPTION callrate 0.90",
  "OPTION callrateAnim 0.80"
)
# file of parameter 
parameter_con <- file("renumf90.par")
writeLines(text = paste(parameter_txt_v, collapse = "\n"),
           con = parameter_con)
close(parameter_con)

#
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
