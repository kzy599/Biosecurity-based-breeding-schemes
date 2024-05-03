#globe parameters

#scheme
sch <- "bs"#ds tds tbs

bsname <- "BSTT"#BSTR BSRR BSRT BSTB BSBT BSBB BSBR BSRB DST TDS TBS
#Gneneration
nGeneration=30

#parents
nFemale = 100 
nMale = 50

#
ChrSize=(2.6 * 10^9) / 44

#population structure
nfamily = 100 #Family number
nProgeny = 200  #the number of progeny per family

#genotype strategy for bio-security selection

#for breeding environment
gsb <- "TOP" #TOP TB RAND

#for commercial environment
gsc<- "TOP" #TOP TB RAND

#genotype strategy for direct selection
gsd <- "TOP"

#the number of candidate family
ncf <- 50

#the number of genotype individual per family
ngi <- 20 # set 10 for easy calculation temporarily

