#setwd("/home/GTDisk1/kangziyi/GxE")
rm(list=ls())
gc()
library(AlphaSimR)
library(data.table)
founderPop = quickHaplo(nInd=1000, nChr=10, segSites=2000, inbred=FALSE)
SP = SimParam$new(founderPop)
SP$restrSegSites(minSnpFreq = 0.05,overlap = TRUE)#设置SNP频率
SP$addTraitAG(100,mean=c(0,0),var = c(1,1),varGxE = c(0.1,0.1),corA = matrix(c(1,0.9,0.9,1),nrow = 2),corGxE = matrix(c(1,0.9,0.9,1),nrow = 2))
SP$setVarE(h2=c(0.41,0.41))#遗传力和偏差
SP$addSnpChip(1250)#55kSNP芯片
SP$setSexes("yes_sys")#按照一个雌一个雄来分配个体的性别
#综上将参数设置好之后进行新的群体的建立，将参数与个体结合
pop_founder = newPop(founderPop, simParam=SP)
cor(pop_founder@gv[,1],pop_founder@gv[,2])

nFemale = 100 #一百个母本
nMale = 50  #50个父本
nCrosses = 100 #产生100个家系
nProgeny = 200  #每个家系200个体
pop <- selectCross(pop_founder,
                   nFemale = nFemale,nMale = nMale,
                   nCrosses = nCrosses,nProgeny = nProgeny,
                   simParam = SP)

getfamgv<- function(pop,popan){
  damid<- unique(pop@mother)
  famgv1 <- c(numeric(length = length(damid)))
  famgv2 <- c(numeric(length = length(damid)))
  for (i in 1:length(damid)) {
    popind<- pop[pop@mother==damid[i]]
    famgv1[i]<- (as.numeric(popan[unique(popind@father)]@gv[,1])/2)+(as.numeric(popan[damid[i]]@gv[,1])/2)
    famgv2[i]<- (as.numeric(popan[unique(popind@father)]@gv[,2])/2)+(as.numeric(popan[damid[i]]@gv[,2])/2)
  }
  return(as.data.table(cbind(famgv1,famgv2)))
}

famgv<- getfamgv(pop = pop,popan = pop_founder)
geresult<- data.table(gs=c(numeric(length =21)),pheno=c(numeric(length =21)))
geresult[1]<- cor(famgv$famgv1,famgv$famgv2)

poppheno<- pop
popgs<- pop
#########
#按照测试环境进行选育，表型和GS对长期GE变化的影响
#gs
pop<- popgs
for (i in 1:20) {
  popan<- pop
  popfam<- selectFam(pop = pop,nFam = 50,trait = 1,use = "gv")
  candidate<- selectWithinFam(pop = popfam,nInd = 6,trait = 1,use = "gv")
  pop <- selectCross(candidate,
                     nFemale = nFemale,nMale = nMale,
                     nCrosses = nCrosses,nProgeny = nProgeny,
                     trait = 1,use = "gv",
                     simParam = SP)
  famgv<- getfamgv(pop = pop,popan = popan)
  geresult$gs[i+1]<- cor(famgv$famgv1,famgv$famgv2)
}
  #pheno
  pop<- poppheno
  for (i in 1:20) {
    popan<- pop
    popfam<- selectFam(pop = pop,nFam = 50,trait = 1,use = "pheno")
    candidate<- selectWithinFam(pop = popfam,nInd = 6,trait = 1,use = "pheno")
    pop <- selectCross(candidate,
                       nFemale = nFemale,nMale = nMale,
                       nCrosses = nCrosses,nProgeny = nProgeny,
                       trait = 1,use = "pheno",
                       simParam = SP)
    famgv<- getfamgv(pop = pop,popan = popan)
    geresult$pheno[i+1]<- cor(famgv$famgv1,famgv$famgv2)
    
  }
  grplot<- cbind(c(0:20),geresult)
  colnames(grplot)[1] <- "Generation"
  dt <- as.data.table(cbind(rep(grplot$Generation,2),c(grplot$gs,grplot$pheno),c(rep("gs",21),rep("pheno",21))))
colnames(dt) <- c("Generation","GxE","method")  
dt$Generation <- as.numeric(dt$Generation)
dt$GxE <- as.numeric(dt$GxE)
picture<- ggplot(data = dt,aes(x=Generation,y=GxE,group=method,linetype=method))+
  geom_point()+
  geom_line()+
  xlab("Generation")+
  ylab("GxE")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text = element_text(family = "STXihei"),
        legend.position = "right",
        legend.background = element_rect(colour = "black"))+
  scale_x_continuous(limits = c(0,20),breaks = seq(0,20,1))
picture                     
