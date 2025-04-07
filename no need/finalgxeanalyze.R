rm(list=ls())
gc()
options(repos = c(CRAN = "https://cran.rstudio.com/"))
library(devtools)
devtools::install_github("tidyverse/ggplot2")
library(ggpubr)
library(purrr)
library(dplyr)
library(data.table)
library(ggpmisc)
library(ggplot2)
library(ggsci)
#提取数据
#在s和p的对应的项目名称处添加 DBST TDBS即可
semodel="s"
if(semodel=="s"){
  gechar<- c("ge2","ge5","ge8")
  for (ge in c(2,5,8)) {
    gelj<- gechar[gechar%flike%ge]
    
    for (genoratio in c(2,5,8)) {
      rm(list=setdiff(ls(),c("genoratio","ge","semodel","gelj","gechar")))
      setwd(paste("/home/GTDisk1/kangziyi/geoutput/",gelj,"/",genoratio,sep = ""))
      programname<- c("BSTB","BSTR","BSTT")#add DBST
      nCycle=10
      #if(ge==2 & genoratio==8){
      # nCycle=3
      #}
      for (p in 1:length(programname)) {
        for(r in 1:nCycle){
            pname<- programname[p]
            out <- fread(input = paste(pname,r,".csv",sep = ""),sep = ",") 
            getWinb<- function(W,idc,Finitial=0,Finb){
              Winb <- W*(1-(idc*(Finitial-Finb)))
              return(Winb)
            }
            getprate = function(x){
              prate = numeric(length = (length(x)-1))
              
              for(i in 2:length(x)){
                prate[i-1] = ((x[i]-x[i-1])/(1-x[i-1]))*100
              }
              prate = rep(mean(prate),(length(x)))
              return(prate)
            }
            idc <- -0
            out$gv<- getWinb(W=out$gv,idc = idc,Finb = out$Inbreeding)
            gg0<- out$gv-out$gv[1]
            #水平值可以这样来分段，
            #如果看不同世代区间的增长率，需要确保基础一致，
            #要用与第0世代的差值作为进展，再用方程同时拟合所有数据再去分段
            gg10<-  out$gv-out$gv[11]
            gg20<-  out$gv-out$gv[21]
            inb <-  out$Inbreeding_plink
            inbp <- out$Inbreeding
            inbrate= getprate(x = inb)
            inbprate = getprate(x = inbp)
            ld = out$LD
            acc <- out$accuracy
            gxxe <- out$gxe
            va1 <- out$va1 #nucleus breeding center
            va2 <- out$va2 #commercial farm environment
            cov12 <- out$cov12
            #gg0 = gg0/sqrt(variance[1])
            #gg10 = gg10/sqrt(variance[11])
            #gg20 = gg20/sqrt(variance[21])
            effsize <- out$Ne
            g2 <- out$g2
            i2 <- out$i2
            covig <- (va1-va2)/4
            nGeneration<- out$generation
            if(r==1){
              gxecounts <- out$gxe[length(nGeneration)]
              mingxecounts <- min(out$gxe)
            }else{
              gxecounts <- c(gxecounts,out$gxe[length(nGeneration)])
              mingxecounts <- c(mingxecounts,min(out$gxe))
            }
            out<- as.data.table(cbind(pname,nGeneration,gg0,gg10,gg20,inbp,inbprate,acc,gxxe,va1,va2,cov12,effsize,g2,i2,covig,ld))
            if(r==1&p==1){
              output<- out
            }else{
              output<- rbind(output,out)
            }
        }
          lcounts<- c(1:length(gxecounts))
          repnum<- lcounts[gxecounts==min(gxecounts)]
          minrepnum<- lcounts[mingxecounts==min(mingxecounts)]
          repdt<- cbind(pname,repnum,minrepnum)
          sdminall <- data.table(pname=pname,sdall=sd(mingxecounts))
          if(p==1){
            repout<- repdt
            sdminall1 <- sdminall
          }else{
            repout<- rbind(repout,repdt)
            sdminall1 <- rbind(sdminall1,sdminall)
          }
            min30th<- fread(input = paste(pname,repnum,".csv",sep = ""),sep = ",")
            min30th[,dttype:=1]
            min30th[,pname:=pname]
            minall <- fread(input = paste(pname,minrepnum,".csv",sep = ""),sep = ",")
            minall[,dttype:=2]
            minall[,pname:=pname]
            minfinal= rbind(min30th,minall)
            if(any(colnames(minfinal)%in%c("numsire", "numdam",  "numid",   "numebv"))){
              minfinal[,numsire:=NULL]
              minfinal[,numdam:=NULL]
              minfinal[,numid:=NULL]
              minfinal[,numebv:=NULL]
            }
            if(p==1){
              minfinal1= minfinal
            }else{
              minfinal1 = rbind(minfinal1,minfinal)
            }
      }
      ##输出目录，output+ge强度+s/p+分型比例
      setwd("/home/GTDisk1/kangziyi/geoutput")
      repout<- as.data.table(repout)
      sdminall1<- as.data.table(sdminall1)
      #其实可以最后再输出，把ge semodel genoratio 都添加到相应的变量里面
      #输出一个文件就足够了
      fwrite(minfinal1,file = paste("minout",ge,semodel,genoratio,"pvalue",".csv",sep = ""))
      fwrite(sdminall1,file = paste("sdminout",ge,semodel,genoratio,"pvalue",".csv",sep = ""))
      fwrite(output,file = paste("output",ge,semodel,genoratio,"pvalue",".csv",sep = ""))
      fwrite(repout,file = paste("repout",ge,semodel,genoratio,"pvalue",".csv",sep = ""))
    }   
  }
}else if(semodel=="p"){
  gechar<- c("ge2","ge5","ge8")
  for(ge in c(2,5,8)){
    gelj<- gechar[gechar%flike%ge]
    genoratio=""
    rm(list=setdiff(ls(),c("ge","semodel","genoratio","gelj","gechar")))
    setwd(paste("/home/GTDisk1/kangziyi/geoutput/",gelj,"/t",sep = ""))
    programname<- c("TBS","TDS","TDBS","DST","DBST")#add TDBS
    nCycle=10
    nrep=3
    
    for (p in 1:length(programname)) {
      for (repc in 1:nrep) {
        for(r in 1:nCycle){
            pname<- programname[p]
            out <- fread(input = paste(pname,r,repc,".csv",sep = ""),sep = ",") 
            getWinb<- function(W,idc,Finitial=0,Finb){
              Winb <- W*(1-(idc*(Finitial-Finb)))
              return(Winb)
            }
            getprate = function(x){
              prate = numeric(length = (length(x)-1))
              
              for(i in 2:length(x)){
                prate[i-1] = ((x[i]-x[i-1])/(1-x[i-1]))*100
              }
              prate = rep(mean(prate),(length(x)))
              return(prate)
            }
            idc <- -0
            out$gv<- getWinb(W=out$gv,idc = idc,Finb = out$Inbreeding)
            gg0<- out$gv-out$gv[1]
            #水平值（绝对值）的变化可以这样来分段，
            #如果看不同世代区间的增长率（相对值），需要确保基础一致，
            #要用与第0世代的差值作为进展，再用方程同时拟合所有数据再去分段
            #线性模型的系数需要相对一个基准来解释，所以分段比较时需要在一个统一基准下
            gg10<- out$gv-out$gv[11]
            gg20<- out$gv-out$gv[21]
            inb <- out$Inbreeding_plink
            inbp<-out$Inbreeding
            inbrate= getprate(x = inb)
            inbprate = getprate(x = inbp)
            ld = out$LD
            acc <- out$accuracy
            gxxe <- out$gxe
            va1 <- out$va1 #nucleus breeding center
            va2 <- out$va2 #commercial farm environment
            cov12 <- out$cov12
            #gg0 = gg0/sqrt(variance[1])
            #gg10 = gg10/sqrt(variance[11])
            #gg20 = gg20/sqrt(variance[21])
            effsize <- out$Ne
            i2 <- out$i2
            g2 <- out$g2
            covig <- (va1-va2)/4
            nGeneration<- out$generation
            if(r==1&repc==1){
              gxecounts <- out$gxe[length(nGeneration)]
              mingxecounts <- min(out$gxe)
            }else{
              gxecounts <- c(gxecounts,out$gxe[length(nGeneration)])
              mingxecounts <- c(mingxecounts,min(out$gxe))
            }
            out<- as.data.table(cbind(pname,nGeneration,gg0,gg10,gg20,inbp,inbprate,acc,gxxe,va1,va2,cov12,effsize,g2,i2,covig,ld))
            if(r==1&p==1&repc==1){
              output<- out
            }else{
              output<- rbind(output,out)
            }
          }
      }
      lcounts<- c(1:length(gxecounts))
      repnum<- lcounts[gxecounts==min(gxecounts)]
      minrepnum<- lcounts[mingxecounts==min(mingxecounts)]
      repdt<- cbind(pname,repnum,minrepnum)
      sdminall <- data.table(pname=pname,sdall=sd(mingxecounts))
      if(p==1){
        repout<- repdt
        sdminall1 <- sdminall
      }else{
        repout<- rbind(repout,repdt)
        sdminall1 <- rbind(sdminall1,sdminall)
      }
      if(repnum<=10){
        rrpc=1
        repnum= repnum
      }else if(repnum<=20){
        rrpc=2
        repnum= repnum-10
      }else{
        rrpc=3
        repnum = repnum -20
      }
      min30th<- fread(input = paste(pname,repnum,rrpc,".csv",sep = ""),sep = ",")
      min30th[,dttype:=1]
      min30th[,pname:=pname]
      if(minrepnum<=10){
        rrpc=1
        minrepnum = minrepnum
      }else if(minrepnum<=20){
        rrpc=2
        minrepnum = minrepnum - 10
      }else{
        rrpc=3
        minrepnum = minrepnum -20
      }
      minall <- fread(input = paste(pname,minrepnum,rrpc,".csv",sep = ""),sep = ",")
      minall[,dttype:=2]
      minall[,pname:=pname]
      minfinal= rbind(min30th,minall)
      if(any(colnames(minfinal)%in%c("numsire", "numdam",  "numid",   "numebv"))){
        minfinal[,numsire:=NULL]
        minfinal[,numdam:=NULL]
        minfinal[,numid:=NULL]
        minfinal[,numebv:=NULL]
      }
      if(p==1){
        minfinal1= minfinal
      }else{
        minfinal1 = rbind(minfinal1,minfinal)
      }
      
    }
    ##输出目录，output+ge强度+s/p+分型比例
    setwd("/home/GTDisk1/kangziyi/geoutput")
    repout<- as.data.table(repout)
    fwrite(minfinal1,file = paste("minout",ge,semodel,genoratio,"pvalue",".csv",sep = ""))
    fwrite(sdminall1,file = paste("sdminout",ge,semodel,genoratio,"pvalue",".csv",sep = ""))
    fwrite(output,file = paste("output",ge,semodel,genoratio,"pvalue",".csv",sep = ""))
    fwrite(repout,file = paste("repout",ge,semodel,genoratio,"pvalue",".csv",sep = ""))
    }
}


rm(list=ls())
gc()
for(semodel in c("s",'p')){
  for(ge in c(2,5,8)){
    if(semodel=="s"){
      setwd("/home/GTDisk1/kangziyi/geoutput")
      outgs2<- fread(paste("output",ge,semodel,2,"pvalue",".csv",sep = ""))
      outgs2[,genoratio:="20%"]
      outgs5<- fread(paste("output",ge,semodel,5,"pvalue",".csv",sep = ""))
      outgs5[,genoratio:="50%"]
      outgs8<- fread(paste("output",ge,semodel,8,"pvalue",".csv",sep = ""))
      outgs8[,genoratio:="80%"]
      #outgs8<- fread(paste("output",ge,semodel,8,"pvalue",".csv",sep = ""))
      #outgs8[,genoratio:="80%"]
      output<- rbind(outgs2,outgs5,outgs8)
      output[,gev:=ge] 
      if(ge==2){
        output1<- output
      }else{
        output1 <-rbind(output1,output)
      }
      
    }else{
      setwd("/home/GTDisk1/kangziyi/geoutput")
      outgs2<- fread(paste("output",ge,semodel,"pvalue",".csv",sep = ""))
      outgs2[!pname%in%c("DST","DBST"),genoratio:="pblup"]
      outgs2[pname%in%c("DST","DBST"),genoratio:="ssgblup"]
      #outgs8<- fread(paste("output",ge,semodel,8,"pvalue",".csv",sep = ""))
      #outgs8[,genoratio:="80%"]
      output<-outgs2
      output[,gev:=ge] 
      if(ge==2){
        outputp1<- output
      }else{
        outputp1 <-rbind(outputp1,output)
      }
    }
  }
}
outputfinal<- rbind(output1,outputp1)

# mymodel<- aov(gg0~gxxe,data= outputfinal[nGeneration==30&pname%in%c("BSTT","BSTB","BSTR")&gev==8,])
# 
# cor(outputfinal[nGeneration==30&pname%in%c("BSTT","BSTB","BSTR")&gev==8,gxxe],
#     outputfinal[nGeneration==30&pname%in%c("BSTT","BSTB","BSTR")&gev==8,gg0])
# meangxe = c()
# for(g in c("20%","50%","80%")){
#   meangxe= c(meangxe,mean(zt1[pname%in%c("BSTB")&gev==5&genoratio==g,gxxe]))
# }
# 
# cor(meangxe,
#     zt1[nGeneration==30&pname%in%c("BSTB")&gev==5,gg0])
# 
# cor(zt1[nGeneration==30&pname%in%c("BSTT","BSTB","BSTR")&gev==2,gxxe],
#     zt1[nGeneration==30&pname%in%c("BSTT","BSTB","BSTR")&gev==2,gg0])
# 
# summary(mymodel)
# lm(log(i2) ~ nGeneration ,data = outputfinal[pname%in%c("BSTT","BSTB","BSTR")&gev==8,])
# lm(log(i2) ~ generation ,data = outputfinal[genoratio=="80%"&gev=="GxE0.2",])
# 
# t.test(x = outputfinal[nGeneration==30&genoratio=="80%",gg0],y=outputfinal[nGeneration==30&genoratio=="50%",gg0])

# for(p in c("BSTR","BSTT")){
#   
#   pvalue = c()
#   
#   for(ge in c(2,5,8)){
#     for(geno in c("20%","50%","80%")){
#       chk = t.test(x = outputfinal[nGeneration==30&gev==ge&genoratio==geno&pname==p,gg0],y=outputfinal[nGeneration==30&gev==ge&genoratio==geno&pname=="BSTB",gg0])
#       
#       pvalue = c(pvalue,chk$p.value)
#     }
#   } 
#   
#   if(p == "BSTR") tr = pvalue else tt = pvalue
# }
# pvalue = data.table(tr = tr, tt = tt)
# pvalue
GEI = 8
GEI = c(2,5,8)
summary(aov(inbprate~pname+genoratio+gev+nGeneration,data =outputfinal[nGeneration%in%c(1:30)&pname%in%c("BSTR","BSTT","BSTB"),] ))

t.test(x = outputfinal[gev %in% GEI&nGeneration==30&pname%in%c("BSTR","BSTT","BSTB"),gg0],y=outputfinal[gev %in% GEI&nGeneration==30&pname=="DST",gg0])
t.test(x = outputfinal[gev %in% GEI&nGeneration==30&pname%in%c("TBS"),gg0],y=outputfinal[gev %in% GEI&nGeneration==30&pname=="TDS",gg0])

t.test(x = outputfinal[gev %in% GEI&nGeneration==30&pname%in%c("BSTR","BSTT","BSTB")&genoratio=="50%",inbprate],y=outputfinal[gev %in% GEI&nGeneration==30&pname=="DST",inbprate])
t.test(x = outputfinal[gev %in% GEI&nGeneration==30&pname%in%c("TBS"),inbprate],y=outputfinal[gev %in% GEI&nGeneration==30&pname=="TDS",inbprate])

t.test(x = outputfinal[gev %in% GEI&nGeneration %in% c(1:30)&pname%in%c("BSTR","BSTT","BSTB")&genoratio=="50%",inbp],y=outputfinal[gev %in% GEI&nGeneration %in% c(1:30)&pname=="DST",inbp])
t.test(x = outputfinal[gev %in% GEI&nGeneration %in% c(1:30)&pname%in%c("TBS"),inbp],y=outputfinal[gev %in% GEI&nGeneration %in% c(1:30)&pname=="TDS",inbp])

t.test(x = outputfinal[gev %in% GEI&nGeneration==30&pname%in%c("BSTR","BSTT","BSTB")&genoratio=="80%",va1],y=outputfinal[gev %in% GEI&nGeneration==30&pname=="DST",va1])
t.test(x = outputfinal[gev %in% GEI&nGeneration==30&pname%in%c("TBS"),va1],y=outputfinal[gev %in% GEI&nGeneration==30&pname=="TDS",va1])

t.test(x = outputfinal[gev %in% GEI&nGeneration==30&pname%in%c("BSTR","BSTT","BSTB")&genoratio=="80%",va2],y=outputfinal[gev %in% GEI&nGeneration==30&pname=="DST",va2])
t.test(x = outputfinal[gev %in% GEI&nGeneration==30&pname%in%c("TBS"),va2],y=outputfinal[gev %in% GEI&nGeneration==30&pname=="TDS",va2])

# dt_test_bs_nbs = outputfinal[nGeneration==30&pname%in%c("BSTR","BSTT","BSTB","DST","TBS","TDS"),]
# dt_test_bs_nbs$gev = factor(dt_test_bs_nbs$gev,levels = c("8","5","2"))
# 
# aov_model_tbs_gev = aov(inbprate ~ gev,data =dt_test_bs_nbs[nGeneration==30&pname%in%c("TBS"),] )
# 
# aov_model_tbs_gev = aov(inbprate ~ gev,data =dt_test_bs_nbs[nGeneration==30&pname%in%c("BSTR","BSTT","BSTB"),] )
# 
# summary(aov_model_tbs_gev)
# 
# tukey_results <- TukeyHSD(aov_model_tbs_gev,"gev",conf.level = 0.95)
# tukey_results
# tukey_df <- broom::tidy(tukey_results)
# # Create the plot using ggplot2
# tukey_plot <- ggplot(tukey_df, aes(x = contrast, y = estimate, ymin = conf.low, ymax = conf.high)) +
#   geom_pointrange() +  # Add point ranges for estimates and confidence intervals
#   geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +  # Add a reference line at y = 0
#   labs(title = "Tukey HSD Test Results", x = "Pairwise Comparisons", y = "Difference in Means") +
#   theme_minimal() +  # Use a minimal theme
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability
# 
# # Print the plot
# print(tukey_plot)
# 
# ggsave("FirstRvision_Figure3.pdf", tukey_plot , width = 15, height = 5, dpi = 300)

# chk = summary(lm(log(gg0)~nGeneration+pname,data = outputfinal[nGeneration%in%c(1:30)&pname%in%c("BSTR","BSTT","BSTB"),]))
# b = chk$coefficients[2,1] + chk$coefficients[3,1] 
# exp(b)
gen = c(1:30)
proname = c("BSTR","BSTT","BSTB","TBS","TDS","DST")
chk = summary(lm(log(1-inbp)~nGeneration+pname+as.factor(gev)+genoratio,data = outputfinal[nGeneration%in%gen&pname%in%proname,]))
b = chk$coefficients[2,1]
b = chk$coefficients[2,1] + chk$coefficients[3,1] + chk$coefficients[10,1]
(1-exp(b))*100

library(ggplot2)
library(dplyr)
library(multcompView)
summary(aov(gg0~pname+genoratio+gev,data =outputfinal[nGeneration==30&pname%in%c("BSTR","BSTT","BSTB"),] ))
dt_test = outputfinal[nGeneration==30&pname%in%c("BSTR","BSTT","BSTB"),]
# colnames(dt_test)[(colnames(dt_test) == "pname")] = ""
dt_test[pname=="BSTR",pname:="RAN"]
dt_test[pname=="BSTB",pname:="T&B"]
dt_test[pname=="BSTT",pname:="TOP"]
dt_test[genoratio=="20%",genoratio:="20"]
dt_test[genoratio=="50%",genoratio:="50"]
dt_test[genoratio=="80%",genoratio:="80"]
dt_test$pname = factor(dt_test$pname,levels = c("T&B","TOP","RAN"))
dt_test$genoratio = factor(dt_test$genoratio,levels = c("80","50","20"))
dt_test$gev = factor(dt_test$gev,levels = c("8","5","2"))
anova_model <- aov(gg0~pname+genoratio+gev,data =dt_test )
anova_model_inbrate <- aov(inbprate~pname+genoratio+gev,data =dt_test )
anova_model_ld <- aov(ld~pname+genoratio+gev,data =dt_test )
summary(anova_model)
summary(anova_model_inbrate)
summary(anova_model_ld)
tukey_results <- TukeyHSD(anova_model,"pname",conf.level = 0.95)
names(tukey_results) = "selective breeding methods for testing group (TG)"
pdf(file = "chk.pdf",width = 15,height = 5)
par(mar = c(5, 8, 4, 2) + 0.1)
plot(tukey_results,las = 1 ,col = "brown")
dev.off()

tukey_results <- TukeyHSD(anova_model,"genoratio",conf.level = 0.95)
names(tukey_results) = "The number of genotyped individuals for selection group (SG)"
pdf(file = "chk2.pdf",width = 15,height = 5)
par(mar = c(5, 8, 4, 2) + 0.1)
plot(tukey_results,las = 1 ,col = "brown")
dev.off()


tukey_results = TukeyHSD(anova_model_inbrate,c("pname","genoratio"),conf.level = 0.95)
tukey_df <- broom::tidy(tukey_results)

# Create the plot using ggplot2
tukey_plot <- ggplot(tukey_df, aes(x = contrast, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_pointrange() +  # Add point ranges for estimates and confidence intervals
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +  # Add a reference line at y = 0
  labs(title = "Tukey HSD Test Results", x = "Pairwise Comparisons", y = "Difference in Means") +
  theme_minimal() +  # Use a minimal theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
        axis.line = element_line(color = "black"),  # Ensure axis lines are visible
        axis.title.y = element_text(color = "black")) +  # Ensure y-axis title is visible
  scale_y_continuous(limits = c(min(tukey_df$conf.low) * 1.1, max(tukey_df$conf.high) * 1.1))  # Adjust y-axis limits
# Print the plot
print(tukey_plot)

#genetic gain
ggsave("FirstRvision_Figure3.pdf", tukey_plot , width = 15, height = 5, dpi = 300)

#inbreeding rate
ggsave("FirstRvision_Figure5.pdf", tukey_plot , width = 15, height = 5, dpi = 300)


summary(lm(log(gg0)~pname+genoratio+factor(gev)+nGeneration,data =outputfinal[nGeneration%in%c(1:30)&pname%in%c("BSTR","BSTT","BSTB"),] ))
summary(aov(gg0~pname+genoratio+gev+nGeneration,data =outputfinal[nGeneration%in%c(1:30)&pname%in%c("BSTR","BSTT","BSTB"),] ))

getzt<- function(output){
  pgratio <- unique(output$genoratio)
  for (p in 1:length(pgratio)) {
    pp<- unique(output[genoratio==pgratio[p],pname])
    for(g in c(0:30)){
      for(i in 1:length(pp)){
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,gg0se:=sd(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$gg0)]
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$gg0 <- mean(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$gg0)
        
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,gg10se:=sd(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$gg10)]
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$gg10 <- mean(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$gg10)
        
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,gg20se:=sd(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$gg20)]
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$gg20 <- mean(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$gg20)
        
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,inbpse:=sd(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$inbp)]
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$inbp <- mean(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$inbp)
        
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,inbpratese:=sd(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$inbprate)]
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$inbprate <- mean(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$inbprate)
        
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,gxese:=sd(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$gxxe)]
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$gxxe <- mean(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$gxxe)
        
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,Nese:=sd(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$effsize)]
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$effsize <- mean(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$effsize)
        
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,va1se:=sd(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$va1)]
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$va1 <- mean(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$va1)
        
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,va2se:=sd(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$va2)]
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$va2 <- mean(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$va2)
        
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,cov12se:=sd(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$cov12)]
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$cov12 <- mean(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$cov12)
        
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,accse:=sd(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$acc)]
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$acc <- mean(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$acc)
        
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,ipgse:=sd(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$i2+output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$g2)]
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,ipg := mean(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$i2+output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$g2)]
        
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,idgse:=sd(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$i2/output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$g2)]
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,idg := mean(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$i2/output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$g2)]
        
        
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,g2se:=sd(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$g2)]
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$g2 <- mean(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$g2)
        
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,i2se:=sd(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$i2)]
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$i2 <- mean(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$i2)
        
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,covigse:=sd(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$covig)]
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$covig <- mean(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$covig)
        
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,ldse:=sd(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$ld)]
        output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$ld <- mean(output[pname==pp[i]&genoratio==pgratio[p]&nGeneration==g,]$ld)
      
        
      }
      
    }
  }
  zt<- unique(output)
  return(zt)
}

#分析均值进展以及近交Ne等
for (ge in c(2,5,8)) {
  if(ge==2){
    zt1<- getzt(outputfinal[gev==ge,.(pname,nGeneration,genoratio,gev,gg0,gg10,gg20,inbp,inbprate,acc,gxxe,va1,va2,cov12,effsize,g2,i2,covig,ld)])
  }else{
    zt1<- rbind(zt1,getzt(outputfinal[gev==ge,.(pname,nGeneration,genoratio,gev,gg0,gg10,gg20,inbp,inbprate,acc,gxxe,va1,va2,cov12,effsize,g2,i2,covig,ld)]))
  }
}

zt<- zt1[genoratio%in%c("80%","pblup","ssgblup")&nGeneration==30,]
zt<- zt1[gev==8&nGeneration==30,]
zt<- zt1[gev==8&nGeneration==30&pname%in%c("TBS","TDS"),]
zt<- zt1[gev==8&nGeneration==30&pname%in%c("BSTT","BSTB","BSTR","DST"),]
zt<- zt1[gev==8&nGeneration==30&pname%in%c("BSTT","BSTB","BSTR","DST","TBS","TDS"),]
zt<- zt1[nGeneration==30&pname%in%c("BSTT","BSTB","BSTR","DST","TBS","TDS"),]
mean(zt$gxese/zt$gxxe)
zt<- zt1[nGeneration==30,]
##作差异柱状图
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent',color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#去网格线
          axis.line = element_line(colour = "black"),
          #axis.title.x = element_blank(),#去x轴标签
          axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
          axis.title.x=element_text(face = "bold",size = 14),
          axis.text = element_text(face = "bold",size = 12),#坐标轴刻度标签加粗
          axis.ticks = element_line(color='black'),
          # axis.ticks.margin = unit(0.8,"lines"),
          #legend.title=element_blank(),
          #legend.position=c(0.5, 0.95),#图例在绘图区域的位置
          #legend.position="none",
          legend.position="right",
          #legend.direction = "horizontal",
          legend.direction = "vertical",
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=8)),
          legend.background = element_rect( linetype="solid",colour ="black")
    )
}
###
library(grid)
correzt = function(...){
zt$gev=as.character(zt$gev)
zt$gev<- paste("GxE0.",zt$gev,sep = "")
zt[genoratio=="ssgblup",genoratio:="50%"]
zt[genoratio=="pblup",genoratio:="PBLUP"]
zt[genoratio!="PBLUP",method:="Genomic selection"]
zt[genoratio=="PBLUP",method:="traditional selection"]
zt[pname=="BSTB",pname:="T&B"]
zt[pname=="BSTR",pname:="RAND"]
zt[pname=="BSTT",pname:="TOP"]
zt[pname=="DST",pname:="NBS"]
zt[pname=="TDS",pname:="NBS"]
zt[pname=="TBS",pname:="BS"]
return(zt)
}
# zt = correzt()
# 
# 
# zt[pname=="TOP",pname:="BS"]
# zt[pname%in%c("BS","NBS")&genoratio=="PBLUP",genoratio:="non-genotyping"]
# P<- ggplot(data = zt[pname%in%c("BS","NBS"),],aes(x = pname, y = gg0, group= genoratio,fill = genoratio))+
#   geom_bar(stat = "identity",position = "dodge")+
#   ylab("Genetic gain")+
#   #xlab("traditional selection")+
#   geom_errorbar(aes(ymax = gg0+gg0se, ymin = gg0-gg0se),
#                 position = position_dodge(0.9), width = 0.15)+##不要下标的话把ymin=mean-se去掉
#   #ylim(0,0.3)+
#   #scale_fill_brewer(palette = "Set1")+
#   theme_zg()+facet_wrap(~ gev)+labs(fill = "genotyping portion\nwithin each family\nof selection population") + theme(plot.title = element_text(hjust = 0.5))
# #+facet_grid(vars(gev),vars(genoratio))
# #+facet_wrap(~ pname)
# 
# P
# ggsave("bsnbs.png", P+xlab("Breeding system") , width = 15, height = 5, dpi = 300)
# 
# 
# #if you want BS to NBS for all genotyping strategy runing the following code
# #just for luan's PPT
# zt[pname=="T&B",met:="Top & Bottom"]
# zt[pname=="TOP",met:="TOP"]
# zt[pname=="RAND",met:="RAND"]
# 
# for( i in c(1:3)){
#   needt=zt[pname=="NBS",]
#   needt2 = zt[pname=="BS",]
#   if(i==1){
#     ordt= zt[!pname%in%c("BS","NBS")]
#   }
#   if(i==1){
#     ordt= rbind(ordt,needt[,met:="Top & Bottom"],needt2[,met:="Top & Bottom"])
#   }else if(i==2){
#     ordt= rbind(ordt,needt[,met:="TOP"],needt2[,met:="TOP"])
#   }else{
#     ordt= rbind(ordt,needt[,met:="RAND"],needt2[,met:="RAND"])
#   }
# }
# 
# zt= ordt
# 
# #if you only want to compared the geneic gain of TOP and NBS,
# #just run the one of the three code
# #else run all the following code
# zt[pname=="T&B",pname:="BS"]
# zt[pname=="TOP",pname:="BS"]
# zt[pname=="RAND",pname:="BS"]
# 
# #run the following code when you wan to compared NBS and BS
# zt[pname%in%c("BS","NBS")&genoratio=="PBLUP",genoratio:="0%"]
# zt[pname%in%c("BS","NBS")&genoratio=="0%",genoratio:="0"]
# zt[pname%in%c("BS","NBS")&genoratio=="20%",genoratio:="20"]
# zt[pname%in%c("BS","NBS")&genoratio=="50%",genoratio:="50"]
# zt[pname%in%c("BS","NBS")&genoratio=="80%",genoratio:="80"]
# zt[pname%in%c("BS","NBS")&gev=="GxE0.2",gev:="GEI 0.2"]
# zt[pname%in%c("BS","NBS")&gev=="GxE0.5",gev:="GEI 0.5"]
# zt[pname%in%c("BS","NBS")&gev=="GxE0.8",gev:="GEI 0.8"]
# 
# #if you want to compared one of the three strategy with the NBS
# #changing the facet(~met+gev) to facet(~gev)
# P<- ggplot(data = zt[pname%in%c("BS","NBS"),],aes(x = pname, y = gg0, group= genoratio,fill = genoratio))+
#   geom_bar(stat = "identity",position = "dodge")+
#   ylab("Genetic gain")+
#   #xlab("traditional selection")+
#   geom_errorbar(aes(ymax = gg0+gg0se, ymin = gg0-gg0se),
#                 position = position_dodge(0.9), width = 0.15)+##不要下标的话把ymin=mean-se去掉
#   #ylim(0,0.3)+
#   #scale_fill_brewer(palette = "Set1")+
#   theme_zg()+facet_wrap(~ gev)+labs(fill = "保种群分型比例\n（%）") + theme(plot.title = element_text(hjust = 0.5))
# #+facet_grid(vars(gev),vars(genoratio))
# #+facet_wrap(~ pname)
# 
# P
# ggsave("bsnbsall.png", P+xlab("Breeding system") , width = 15, height = 5, dpi = 300)
#==================作图真正的开始=================================================================================
#========================================================================================================================
#the loss of genetic gain
#genomic 

# gg0_tb = zt1[nGeneration == 30&pname=="BSTB"&gev%in%c(2,5,8)&genoratio%in%c("20%","50%","80%"),gg0] %>% round(2)
# gg0_top = zt1[nGeneration == 30&pname=="BSTT"&gev%in%c(2,5,8)&genoratio%in%c("20%","50%","80%"),gg0] %>% round(2)
# gg0_ran = zt1[nGeneration == 30&pname=="BSTR"&gev%in%c(2,5,8)&genoratio%in%c("20%","50%","80%"),gg0] %>% round(2)

gg0_tb = c(65.60, 74.86, 79.70, 73.91, 80.82, 82.96, 79.89, 83.52, 83.86)
gg0_top = c(64.59, 73.38, 77.27, 71.93, 78.69, 80.59, 78.69, 81.28, 81.87)
gg0_ran = c(64.57, 73.09, 76.48, 72.23, 78.55, 80.04, 78.04, 81.74, 81.49)
gg0_bs = c(gg0_tb,gg0_top,gg0_ran)

GEI = c(rep("GEI 0.2",3),rep("GEI 0.5",3),rep("GEI 0.8",3))
geno = rep(c("20","50","80"),3)
cnum = length(GEI)
bsname = c(rep("T&B",cnum),rep("TOP",cnum),rep("RAN",cnum))

zt_bsg = data.table(pname = bsname, gev = rep(GEI,3),genoratio = rep(geno,3),method = rep("Genomic selection",3*cnum),gg0 = gg0_bs)

#pedigree
# BS_Pedigree = zt1[nGeneration == 30&pname=="TBS"&gev%in%c(2,5,8),gg0] %>% round(2)
BS_Pedigree = c(45.36, 57.16, 68.45)
zt_bsp = data.table(pname = rep("PED",3), gev = c("GEI 0.2","GEI 0.5","GEI 0.8"),genoratio = rep("0",3),method = rep("Pedigree-based selection",3),gg0 = BS_Pedigree)

zt = rbind(zt_bsp,zt_bsg)

#covert to percentage
# NBS_Pedigree = zt1[nGeneration == 30&pname=="TDS"&gev%in%c(2,5,8),gg0] %>% round(2)
# NBS_Genomic = zt1[nGeneration == 30&pname=="DST"&gev%in%c(2,5,8),gg0] %>% round(2)
NBS_Pedigree = c(74.27, 76.03,75.52)
NBS_Genomic = c(85.15, 86.59, 85.69)
GEI = unique(zt$gev)
for(ii in c(1:3)){
  zt[gev==GEI[ii]&method=="Genomic selection",gg0:=round(((gg0-NBS_Genomic[ii])/NBS_Genomic[ii])*100,2)]
  zt[gev==GEI[ii]&method=="Pedigree-based selection",gg0:=round(((gg0-NBS_Pedigree[ii])/NBS_Pedigree[ii])*100,2)]
}

zt
#
GEI = c("GEI 0.2","GEI 0.5","GEI 0.8")
geno = c("50","80")
min(zt[gev %in% GEI & genoratio %in% geno & method == "Genomic selection",gg0])
max(zt[gev %in% GEI & genoratio %in% geno & method == "Genomic selection",gg0])
mean(zt[gev%in% GEI & genoratio %in% geno & method == "Genomic selection",gg0]) %>% round(2)
zt[gev== GEI & method == "Pedigree-based selection",gg0]
zt[gev== "GEI 0.8" & genoratio == "20" & pname == "T&B" & method == "Genomic selection",gg0]
#
#
zt$pname = factor(zt$pname,levels = c("PED","TOP","RAN","T&B"))
zt$width = ifelse(zt$pname == "PED",0.9,0.9)
P<- ggplot(data = zt,aes(x = pname, y = gg0, group=genoratio,fill = genoratio,width = width))+
  geom_bar(stat = "identity",position = "dodge")+#scale_fill_discrete(breaks = c(0,20,50,80))+
  # geom_text(data = subset(zt, genoratio %in% c("0", "20","50","80")), 
  #           aes(label = gg0), 
  #           position = position_stack(vjust = 1),
  #           size = 2,
  #           color = "black") +
  ylab("Loss of genetic gain")+
  theme_zg()+facet_wrap(~ gev)+ theme(plot.title = element_text(hjust = 0.5))+
  labs(fill = "Genotyping proprotions\n within candidate familes\n of breeding popupation")+
  xlab("")+labs(fill = "Number of\nGenotyped\nSG individuals\n(per family)")+
  # theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())
  theme(plot.title = element_text(hjust = 0.5))
P
#the loss of genetic gain
ggsave("FirstRvision_Figure2.pdf", P , width = 15, height = 5, dpi = 300)

##Selective genotyping method within each candidate family of TG
TBgg = zt_bsg[pname == "T&B" & gev %in% c("GEI 0.2","GEI 0.5","GEI 0.8") & genoratio%in%c("20","50","80"),.(gg0,pname,gev,genoratio)]
TOPgg = zt_bsg[pname == "TOP" & gev %in% c("GEI 0.2","GEI 0.5","GEI 0.8") & genoratio%in%c("20","50","80"),.(gg0,pname,gev,genoratio)]
RANgg = zt_bsg[pname == "RAN" & gev %in% c("GEI 0.2","GEI 0.5","GEI 0.8") & genoratio%in%c("20","50","80"),.(gg0,pname,gev,genoratio)]
chk = round(((TBgg$gg0 - TOPgg$gg0)/TOPgg$gg0)*100,2)
min(chk)
max(chk)
mean(chk,2)
chk = round(((TBgg$gg0 - RANgg$gg0)/TOPgg$gg0)*100,2)
#

##The number of genotyped individuals within each candidate family of SG
gg20 = zt_bsg[genoratio == "20" & gev %in% c("GEI 0.2","GEI 0.5","GEI 0.8") & pname%in%c("T&B","TOP","RAN"),.(gg0,pname,gev,genoratio)]
gg50 = zt_bsg[genoratio == "50" & gev %in% c("GEI 0.2","GEI 0.5","GEI 0.8") & pname%in%c("T&B","TOP","RAN"),.(gg0,pname,gev,genoratio)]
gg80 = zt_bsg[genoratio == "80" & gev %in% c("GEI 0.2","GEI 0.5","GEI 0.8") & pname%in%c("T&B","TOP","RAN"),.(gg0,pname,gev,genoratio)]
chk = round(((gg50$gg0 - gg20$gg0)/gg20$gg0)*100,2)
min(chk)
max(chk)
mean(chk,2)
chk = round(((gg80$gg0 - gg50$gg0)/gg50$gg0)*100,2)

GEI = "GEI 0.8"
gg20_GEI = zt_bsg[genoratio == "20" & gev == GEI & pname%in%c("T&B","TOP","RAN"),.(gg0,pname,gev,genoratio)]
gg50_GEI = zt_bsg[genoratio == "50" & gev == GEI & pname%in%c("T&B","TOP","RAN"),.(gg0,pname,gev,genoratio)]
gg80_GEI = zt_bsg[genoratio == "80" & gev == GEI & pname%in%c("T&B","TOP","RAN"),.(gg0,pname,gev,genoratio)]
chk = round(((gg50_GEI$gg0 - gg20_GEI$gg0)/gg20_GEI$gg0)*100,2)
min(chk)
max(chk)
mean(chk,2)
chk = round(((gg80_GEI$gg0 - gg50_GEI$gg0)/gg50_GEI$gg0)*100,2)
#

#Genomic selection vs pedigree-based selection

#BS
GEI = c("GEI 0.2","GEI 0.5","GEI 0.8")
chk = c()
for(i in c(1:3)){
  chk = c(chk,round(((zt_bsg[gev == GEI[i],gg0] - BS_Pedigree[i])/BS_Pedigree[i])*100,2))
}
min(chk)
max(chk)
mean(chk,2)

#NBS
chk = round(((NBS_Genomic - NBS_Pedigree)/NBS_Pedigree)*100,2)
min(chk)
max(chk)
mean(chk,2)
#

##The GEI intensity changes over generations
cvalue = c(numeric(3))
BDschemes = c("TBS")
BDschemes = c("BSTT","BSTB","BSTR")
gen = 30
for(ii in c(2,5,8)){
  MGE = (zt1[gev==ii&nGeneration%in%c(gen:30)&pname%in%BDschemes,.(pname,gxxe,gxese)]$gxxe %>% round(2))
  SDE = (zt1[gev==ii&nGeneration%in%c(gen:30)&pname%in%BDschemes,.(pname,gxxe,gxese)]$gxese %>% round(2))
  cvalue[which(ii == c(2,5,8))] = round(((SDE/MGE)*100),2) %>% mean(2)
}
cvalue
#

# Rate of inbreeding
BDschemes = c("BSTT","BSTB","BSTR","TBS","DST","TDS")
chk = zt1[pname%in%BDschemes&nGeneration==30,.(pname,gev,genoratio,round(inbprate,2),round(inbpratese,2))]
chk$pname = factor(chk$pname,levels =c("BSTB","BSTT","BSTR","DST","TBS","TDS") )
setorder(chk,pname)
chk

chk[genoratio=="50%"&gev%in%c(2,5,8)&pname%in%c("BSTT","BSTB","BSTR"),V4] %>% mean(2)
chk[genoratio%in%c("20%","50%","80%")&gev%in%c(2,5,8)&pname=="BSTT",V4] %>% mean(2)

chk[pname%in%c("BSTT","BSTB","BSTR"),V4] %>% mean(2)
chk[pname%in%c("DST"),V4] %>% mean(2)

chk[pname%in%c("TBS"),V4] %>% mean(2)

chk[pname%in%c("TDS"),V4] %>% mean(2)
#

zt<- zt1[nGeneration==30&pname%in%c("BSTT","BSTB","BSTR","DST","TBS","TDS"),]
zt = correzt()
##just compared samping methopd
zt[genoratio=="PBLUP",genoratio:="0"]
zt[genoratio=="20%",genoratio:="20"]
zt[genoratio=="50%",genoratio:="50"]
zt[genoratio=="80%",genoratio:="80"]
zt[gev=="GxE0.2",gev:="GEI 0.2"]
zt[gev=="GxE0.5",gev:="GEI 0.5"]
zt[gev=="GxE0.8",gev:="GEI 0.8"]
zt$pname=factor(zt$pname,levels = c("BS","RAND","TOP","T&B","NBS"))
zt = zt[method=="Genomic selection"&pname%in%c("RAND","TOP","T&B"),]
zt[pname=="RAND",pname:="RAN"]
zt$pname=factor(zt$pname,levels = c("RAN","TOP","T&B"))
P<- ggplot(data = zt,aes(x = genoratio, y = gg0, group=pname,fill = pname))+
  geom_bar(stat = "identity",position = "dodge")+
  ylab("Genetic gain")+
  geom_errorbar(aes(ymax = gg0+gg0se, ymin = gg0-gg0se),
                position = position_dodge(0.9), width = 0.15)+##不要下标的话把ymin=mean-se去掉
  theme_zg()+facet_wrap(~ gev)+ theme(plot.title = element_text(hjust = 0.5))+#ggtitle("Biosecurity-based breeding schemes with genomic selection") +
  labs(fill = "Genotyping proprotions\n within candidate familes\n of breeding popupation")+
  theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())+xlab("")
P
ggsave("allsampling.png", P , width = 15, height = 5, dpi = 300)
ggsave("Figure3.pdf", P , width = 15, height = 5, dpi = 300)
##

##compared BS to NBS , BS is T&B
#scale_fill_discrete(breaks = c(0,20,50,80))
#no ggtitle
#facet_wrap(~ method+gev)
#if want to compared GS to PED need run the additional code
zt<- zt1[nGeneration==30&pname%in%c("BSTT","BSTB","BSTR","DST","TBS","TDS"),]
zt = correzt()
zt[genoratio=="PBLUP",genoratio:="0"]
zt[genoratio=="20%",genoratio:="20"]
zt[genoratio=="50%",genoratio:="50"]
zt[genoratio=="80%",genoratio:="80"]
zt[gev=="GxE0.2",gev:="GEI 0.2"]
zt[gev=="GxE0.5",gev:="GEI 0.5"]
zt[gev=="GxE0.8",gev:="GEI 0.8"]
zt$pname=factor(zt$pname,levels = c("BS","RAND","TOP","T&B","NBS"))
zt= zt[pname%in%c("BS","T&B","NBS"),]
zt[pname=="T&B",pname:="BS"]
zt[pname=="BS"&method=="traditional selection",genoratio:=0]
zt[pname=="NBS",genoratio:=NA]
zt[method=="traditional selection",method:="Pedigree-based selection"]
P<- ggplot(data = zt,aes(x = pname, y = gg0, group=genoratio,fill = genoratio))+
  geom_bar(stat = "identity",position = "dodge")+scale_fill_discrete(breaks = c(0,20,50,80))+
  ylab("Genetic gain")+
  geom_errorbar(aes(ymax = gg0+gg0se, ymin = gg0-gg0se),
                position = position_dodge(0.9), width = 0.15)+##不要下标的话把ymin=mean-se去掉
  theme_zg()+facet_wrap(~ method+gev)+ theme(plot.title = element_text(hjust = 0.5))+
  labs(fill = "Genotyping proprotions\n within candidate familes\n of breeding popupation")+
  theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())+xlab("")
P
#BS vs NBS
ggsave("Figure2.pdf", P , width = 15, height = 5, dpi = 300)
#additional code, run Figure2 firstly
zt[method=="Pedigree-based selection",method:="Pedigree"]
zt[method=="Genomic selection",method:="Genomic"]
zt$method = factor(zt$method,levels = c("Pedigree","Genomic"))
P<- ggplot(data = zt,aes(x = method, y = gg0, group=genoratio,fill = genoratio))+
  geom_bar(stat = "identity",position = "dodge")+scale_fill_discrete(breaks = c(0,20,50,80))+
  ylab("Genetic gain")+
  geom_errorbar(aes(ymax = gg0+gg0se, ymin = gg0-gg0se),
                position = position_dodge(0.9), width = 0.15)+##不要下标的话把ymin=mean-se去掉
  theme_zg()+facet_wrap(~ pname+gev)+ theme(plot.title = element_text(hjust = 0.5))+
  labs(fill = "Genotyping proprotions\n within candidate familes\n of breeding popupation")+xlab("")+
  theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())
P
#ggsave("GSvsPS.png", P , width = 15, height = 5, dpi = 300)
ggsave("Figure4.pdf", P , width = 15, height = 5, dpi = 300)
#====================================================================
#bs vs nbs 只要pedigree-based作为P1 zt[method== "Pedigree-based selection",]
#ped vs genomic 只要BS作为P2 zt[pname== "BS",]，在P2的theme()中添加legend.direction = "horizontal"
P1
P2
library(ggpubr)
library(ggplot2)
ggarrange(P2+ylim(0,88), P1+ylim(0,88), common.legend = TRUE, legend="top")


#====================================================================
# P1<- ggplot(data = zt[method=="Genomic selection",],aes(x = pname, y = gg0, group= genoratio,fill = genoratio))+
#   geom_bar(stat = "identity",position = "dodge")+
#   ylab("Genetic gain")+
#   #xlab("traditional selection")+
#   geom_errorbar(aes(ymax = gg0+gg0se, ymin = gg0-gg0se),
#                 position = position_dodge(0.9), width = 0.15)+##不要下标的话把ymin=mean-se去掉
#   #ylim(0,0.3)+
#   #scale_fill_brewer(palette = "Set1")+
#   theme_zg()+facet_wrap(~ gev)+ theme(plot.title = element_text(hjust = 0.5))+ggtitle("Genomic selection") +
#   labs(fill = "Genotyping proprotions\n within candidate familes\n of breeding popupation")+xlab("Breeding schemes")+
#   theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())
# #+facet_grid(vars(gev),vars(genoratio))
# #+facet_wrap(~ pname)
# 
# P1
# ggsave("allbs.png", P1 , width = 15, height = 5, dpi = 300)
# 
# 
# P2<- ggplot(data = zt[method=="traditional selection",],aes(x = pname, y = gg0, group= genoratio))+
#   geom_bar(stat = "identity",position = "dodge")+
#   ylab("Genetic gain")+
#   #xlab("traditional selection")+
#   geom_errorbar(aes(ymax = gg0+gg0se, ymin = gg0-gg0se),
#                 position = position_dodge(0.9), width = 0.15)+##不要下标的话把ymin=mean-se去掉
#   #ylim(0,0.3)+
#   #scale_fill_brewer(palette = "Set1")+
#   theme_zg()+facet_wrap(~ gev)+ggtitle("Traditional selection") + theme(plot.title = element_text(hjust = 0.5))
# #+facet_grid(vars(gev),vars(genoratio))
# #+facet_wrap(~ pname)
# 
# P2
# 
# mymodel<- aov(gg0~pname+genoratio,data= zt)
# 
# summary(mymodel)
# t.test(x = outputfinal[nGeneration==30&pname=="BSTT",gg0],y=outputfinal[nGeneration==30&pname=="BSTR",gg0])
# 
# t.test(x = outputfinal[nGeneration==30&pname=="BSTR"&genoratio=="80%",gg0],y=outputfinal[nGeneration==30&pname=="BSTR"&genoratio=="20%",gg0])
# t.test(x = outputfinal[nGeneration==30&pname=="BSTR"&genoratio=="80%",gg0],y=outputfinal[nGeneration==30&pname=="BSTR"&genoratio=="50%",gg0])
# 
# figure <- ggarrange(P1+labs(fill="genotyping portion\nwithin each family") +ylab("Genetic gain")+rremove("xlab")+ylim(0,40) , P2 + rremove("xlab")+ylab("Genetic gain")+ylim(0,40) , # remove axis labels from plots
#                     labels = c("a","b"),
#                     ncol = 1, nrow =2,
#                     widths = c(1,1),
#                     common.legend = TRUE, legend = "right",
#                     align = "v", 
#                     font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
# figure
# ggsave("gxegg.png", figure , width = 10, height = 5, dpi = 300)
# #==================================================================================================================================
# 
# #======================================================================================
# zt<- zt1[nGeneration==30&pname%in%c("BSTT","BSTB","BSTR","DST","TBS","TDS"),]
# zt = correzt()
# zt[genoratio=="PBLUP",genoratio:="0"]
# zt[genoratio=="20%",genoratio:="20"]
# zt[genoratio=="50%",genoratio:="50"]
# zt[genoratio=="80%",genoratio:="80"]
# zt[gev=="GxE0.2",gev:="GEI 0.2"]
# zt[gev=="GxE0.5",gev:="GEI 0.5"]
# zt[gev=="GxE0.8",gev:="GEI 0.8"]
# zt$pname=factor(zt$pname,levels = c("BS","RAND","TOP","T&B","NBS"))
# P1<- ggplot(data = zt,aes(x = pname, y = gg0, group= genoratio,fill = genoratio))+
#   geom_bar(stat = "identity",position = "dodge")+
#   ylab("Genetic gain")+
#   #xlab("traditional selection")+
#   geom_errorbar(aes(ymax = gg0+gg0se, ymin = gg0-gg0se),
#                 position = position_dodge(0.9), width = 0.15)+##不要下标的话把ymin=mean-se去掉
#   #ylim(0,0.3)+
#   #scale_fill_brewer(palette = "Set1")+
#   #theme_zg()+facet_wrap(~ gev)+ theme(plot.title = element_text(hjust = 0.5))+
#   theme_zg()+facet_wrap(~ gev)+ theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())+
#   #labs(fill = "The number of\n genotyped individuals\n within each candidate family\n of breeding population")+
#   xlab("Breeding schemes")
# #+facet_grid(vars(gev),vars(genoratio))
# #+facet_wrap(~ pname)
# 
# P1
# ggsave("ALL.png", P1 , width = 15, height = 5, dpi = 300)
# #====================================================================================================
# 
# 
# zt<- zt1[pname=="TBS",]
# zt<- zt1[genoratio=="80%",]
# zt<- zt1[pname%in%c("BSTT","BSTB","BSTR","TBS"),]
# zt[gev==2&nGeneration==30,]
# zt$gev<- as.character(zt$gev)
# zt$gev<- paste("GxE0.",zt$gev,sep = "")
# zt$pname <- paste(zt$pname,"_",zt$gev,sep = "")
# picture<- ggplot(data = zt,aes(nGeneration,y=gxxe,group=pname,color=gev))+
#   geom_point()+
#   geom_line()+
#   xlab("Generation")+
#   #scale_shape_manual(values = c(3,15,0,NA))+
#   ylab("GxE interaction")+
#   theme_bw()+
#   theme(panel.grid.major = element_line(colour = NA),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank(),
#         text = element_text(family = "STXihei"),
#         legend.position = "right",
#         legend.background = element_rect(colour = "black"))+
#   scale_x_continuous(limits = c(0,30),breaks = seq(0,30,5))+facet_wrap(~ gev+genoratio)+ggtitle("biosecurity schemes") + theme(plot.title = element_text(hjust = 0.5))
# picture


rm(list=ls())
gc()
setwd("/home/GTDisk1/kangziyi/geoutput/")
# semodel="p"
#提取各ge各分型比例下，各方案得到最剧烈一次gxe
#把repout改成minout即可
#再改一下名字
#按pblup和ssgblup分组即可
# outname = "repout"
outname = "minout"
# outname = "sdminout"
#如果是直接输出的一个文件，则直接读取一个文件就行，并不需要分别读取再合并了
for (semodel in c("s","p")) {
  setwd("/home/GTDisk1/kangziyi/geoutput/")
  for (ge in c(2,5,8)) { 
      if(semodel=="s"){
        outgxe2<- fread(paste(outname,ge,semodel,2,"pvalue",".csv",sep = ""),sep = ",")
        outgxe2[,genoratio:="20%"]
        outgxe5<- fread(paste(outname,ge,semodel,5,"pvalue",".csv",sep = ""),sep = ",")
        outgxe5[,genoratio:="50%"]
        outgxe8<- fread(paste(outname,ge,semodel,8,"pvalue",".csv",sep = ""),sep = ",")
        outgxe8[,genoratio:="80%"]
        output<- rbind(outgxe2,outgxe5,outgxe8)
        output[,gev:=ge] 
        if(ge==2){
          gxeput1<- output
        }else{
          gxeput1 <-rbind(gxeput1,output)
        }
      }else{
        outgxe2<- fread(paste(outname,ge,semodel,"pvalue",".csv",sep = ""),sep = ",")
        outgxe2[!pname%in%c("DST","DBST"),genoratio:="pblup"]
        outgxe2[pname%in%c("DST","DBST"),genoratio:="ssgblup"]
        output<- outgxe2
        output[,gev:=ge] 
        if(ge==2){
          gxeputp1<- output
        }else{
          gxeputp1 <-rbind(gxeputp1,output)
        }
      }
    }
  }
gxeputfinal<- rbind(gxeput1,gxeputp1)
#gxeputfinal<- gxeputp1

#如果是minout，则这部分可以省略
# gechar<- c("ge2","ge5","ge8")
# for (ge in c(2,5,8)) {
#   gelj<- gechar[gechar%flike%ge]
#   geno<- c("pblup","ssgblup",2,5,8)
#   #geno<- "pblup"
#   for (g in 1:length(geno)) {
#     if(geno[g]%in%c("pblup","ssgblup")){
#       setwd(paste("/home/GTDisk1/kangziyi/geoutput/",gelj,"/t",sep = ""))
#       programname<- unique(gxeputfinal[gev==ge&genoratio%flike%geno[g],pname])
#       for (p in 1:length(programname)) {
#         rep30min<- gxeputfinal[gev==ge&pname==programname[p],repnum]
#         repallmin <- gxeputfinal[gev==ge&pname==programname[p],minrepnum]
#         if(rep30min<=10){
#           out <- fread(input = paste(programname[p],rep30min,1,".csv",sep = ""),sep = ",")
#           out[,gev:=ge]
#           out[,genorattion:=geno[g]]
#           out[,pname:=programname[p]]
#         }else if(10<rep30min&rep30min<=20){
#           out <- fread(input = paste(programname[p],(rep30min-10),2,".csv",sep = ""),sep = ",")
#           out[,gev:=ge]
#           out[,genorattion:=geno[g]]
#           out[,pname:=programname[p]]
#         }else if(20<rep30min&rep30min<=30){
#           out <- fread(input = paste(programname[p],(rep30min-20),3,".csv",sep = ""),sep = ",")
#           out[,gev:=ge]
#           out[,genorattion:=geno[g]]
#           out[,pname:=programname[p]]
#         }
#         
#         if(repallmin<=10){
#           outall <- fread(input = paste(programname[p],repallmin,1,".csv",sep = ""),sep = ",")
#           outall[,gev:=ge]
#           outall[,genorattion:=geno[g]]
#           outall[,pname:=programname[p]]
#         }else if(10<repallmin&repallmin<=20){
#           outall <- fread(input = paste(programname[p],(repallmin-10),2,".csv",sep = ""),sep = ",")
#           outall[,gev:=ge]
#           outall[,genorattion:=geno[g]]
#           outall[,pname:=programname[p]]
#         }else if(20<repallmin&repallmin<=30){
#           outall <- fread(input = paste(programname[p],(repallmin-20),3,".csv",sep = ""),sep = ",")
#           outall[,gev:=ge]
#           outall[,genorattion:=geno[g]]
#           outall[,pname:=programname[p]]
#         }
#         #out[,dttype:=1]
#         #outall[.dttype:=2]
#         outgxe<- rbind(out,outall)
#         if(ge==2&p==1&g==1){
#           outgxep<- outgxe
#         }else{
#           outgxep<- rbind(outgxep,outgxe)
#         }
#       }
#     }else{
#       setwd(paste("/home/GTDisk1/kangziyi/geoutput/",gelj,"/",geno[g],sep = ""))
#       programname<- unique(gxeputfinal[gev==ge&genoratio%flike%geno[g],pname])
#       for (p in 1:length(programname)) {
#         rep30min<- gxeputfinal[gev==ge&genoratio%flike%geno[g]&pname==programname[p],repnum]
#         repallmin <- gxeputfinal[gev==ge&genoratio%flike%geno[g]&pname==programname[p],minrepnum]
#         
#         out <- fread(input = paste(programname[p],rep30min,".csv",sep = ""),sep = ",")
#         out[,gev:=ge]
#         out[,genorattion:=geno[g]]
#         out[,pname:=programname[p]]
#         
#         outall <- fread(input = paste(programname[p],repallmin,".csv",sep = ""),sep = ",")
#         outall[,gev:=ge]
#         outall[,genorattion:=geno[g]]
#         outall[,pname:=programname[p]]
#         
#         #out[,dttype:=1]
#         #outall[.dttype:=2]
#         
#         outgxe<- rbind(out,outall)
#         if(p==1&ge==2&g==3){
#           outgxe = outgxe[,.(pname,generation,genorattion,gev,gxe)]
#           outgxes<- outgxe
#         }else{
#           outgxe= outgxe[,.(pname,generation,genorattion,gev,gxe)]
#           outgxes<- rbind(outgxes,outgxe)
#         }
#       }
#     }
#   }
# }
# 
# unique(outgxep$pname)
# min(outgxep[pname=="TBS"&gev==8,][1:31,gxe])
# min(outgxep[pname=="TBS"&gev==8,][32:62,gxe])
# min(outgxes[pname=="BSTB"&gev==8&genorattion==2,][1:31,gxe])
# min(outgxes[pname=="BSTB"&gev==8,][32:62,gxe])
# min(outgxes[pname=="BSTR"&gev==8,][1:31,gxe])
# min(outgxes[pname=="BSTR"&gev==8,][32:62,gxe])
# min(outgxes[pname=="BSTT"&gev==8,][1:31,gxe])
# min(outgxes[pname=="BSTT"&gev==8,][32:62,gxe])
#如果添加了dttype则直接改名字
#如果是minout，则需要用
finaldt=gxeputfinal[pname%in%c("BSTT","BSTB","BSTR","TBS"),]
#如果是repout则
#finaldt = outgxep 
#finaldt$pname=paste(finaldt$pname,finaldt$gev,finaldt$dttype,finaldt$genoratio,sep = "_")
dtplot = finaldt

#如果是repout且没有添加dttype则需要下列步骤
# getplotdt_p= function(outgxep,p){
#   for(g in c(2,5,8)){
#     if(g==2){
#       plotdt1<- outgxep[pname==p&gev==g,][1:31,.(gxe,pname,generation,gev,genorattion)]
#     }else{
#       plotdt<- outgxep[pname==p&gev==g,][1:31,.(gxe,pname,generation,gev,genorattion)]
#       plotdt1<- rbind(plotdt1,plotdt)
#     }
#   }
#   plotdt1[,dttype:=1]
#   
#   for(g in c(2,5,8)){
#     if(g==2){
#       plotdt2<- outgxep[pname==p&gev==g,][32:62,.(gxe,pname,generation,gev,genorattion)]
#     }else{
#       plotdt<- outgxep[pname==p&gev==g,][32:62,.(gxe,pname,generation,gev,genorattion)]
#       plotdt2<- rbind(plotdt2,plotdt)
#     }
#   }
#   plotdt2[,dttype:=2]
#   finaldt<- rbind(plotdt1,plotdt2)
#   finaldt$pname=paste(finaldt$pname,finaldt$gev,finaldt$dttype,finaldt$genorattion,sep = "_")
#   return(finaldt)
#   
# }
# getplotdt_s<- function(outgxep,p){
#   for(geno in c(2,5,8)){
#     for(g in c(2,5,8)){
#       if(g==2){
#         plotdt1<- outgxep[pname==p&gev==g&genorattion==geno,][1:31,.(gxe,pname,generation,gev,genorattion)]
#       }else{
#         plotdt<- outgxep[pname==p&gev==g&genorattion==geno,][1:31,.(gxe,pname,generation,gev,genorattion)]
#         plotdt1<- rbind(plotdt1,plotdt)
#       }
#     }
#     plotdt1[,dttype:=1]
#     
#     for(g in c(2,5,8)){
#       if(g==2){
#         plotdt2<- outgxep[pname==p&gev==g&genorattion==geno,][32:62,.(gxe,pname,generation,gev,genorattion)]
#       }else{
#         plotdt<- outgxep[pname==p&gev==g&genorattion==geno,][32:62,.(gxe,pname,generation,gev,genorattion)]
#         plotdt2<- rbind(plotdt2,plotdt)
#       }
#     }
#     plotdt2[,dttype:=2]
#     
#     if(geno==2){
#       finaldt<- rbind(plotdt1,plotdt2)
#       finaldt$pname=paste(finaldt$pname,finaldt$gev,finaldt$dttype,finaldt$genorattion,sep = "_")
#     }else{
#       finaldt1<- rbind(plotdt1,plotdt2)
#       finaldt1$pname=paste(finaldt1$pname,finaldt1$gev,finaldt1$dttype,finaldt1$genorattion,sep = "_")
#       finaldt= rbind(finaldt,finaldt1)
#     }
#   }
#   
#   return(finaldt)
#   
# }
# dtplot<- getplotdt_s(outgxep = outgxes,p = "BSTB")
# dtplot<- getplotdt_p(outgxep = outgxep,p = "TBS")
# dtplot$gev<- as.character(dtplot$gev)


#dtplot1[generation==26,]
#generation,y=gxe,group=pname,color=dttype
#也可以改用三线表格的形式，gxe最大一次降低（取那个世代的se或者重复间30世代最小的一次的标准差），和最后一个世代的降低（重复间具有se的，zt中有）
#需要注意outputfinal是否改变了
#getminse<- function(output){
#pgratio <- unique(output$genoratio)
#for (p in 1:length(pgratio)) {
  #pp<- unique(output[genoratio==pgratio[p],pname])
    #for(i in 1:length(pp)){
      #for(r in 1:10){
        #outp = output[((31*r)-30):(31*r),]
        #outp[pname==pp[i]&genoratio==pgratio[p],mingxe:=min(outp[pname==pp[i]&genoratio==pgratio[p],]$gxxe)]
       #outp[pname==pp[i]&genoratio==pgratio[p],gxe30:=outp[generation==30&pname==pp[i]&genoratio==pgratio[p],]$gxxe]
       # if(r==1&i==1&p==1){
      #    outputneed= outp
     #   }else{
   #       outputneed= rbind(outputneed,outp)
    #    }
  #    }
      
 #   }
    

#}
#zt<- unique(outputneed[.(pname,genoratio,gev,mingxe,gxe30)])
#return(zt)
#}
#for (ge in c(2,5,8)) {
 # if(ge==2){
  #  ztmin<- getminse(outputfinal[gev==ge,])
  #}else{
  #  ztmin<- rbind(ztmin,getmin(outputfinal[gev==ge,]))
  #}
#}
#ztsd= unique(ztmin[pname=="BSTT"&genoratio==2&gev==2,.(pname,mingxe,gxe30)])
#sd(ztsd$mingxe)
#min(ztsd$mingxe)
#sd(ztsd$gxe30)
#min(ztsd$gxe30)

#zt= zt1[pname%in%c("BSTT","BSTB","BSTR","TBS"),]
#dtplot1=dtplot[dttype==1,]
#for (p in c("BSTT","BSTB","BSTR","TBS")) {
 # for(g in c("20%","50%","80%","pblup")){
#    for (v in c(2,5,8)) {
#      dtplot1[pname==p&gev==v&genoratio==g,gxese:=zt[pname==p&gev==v&genoratio==g,gxese]]
#    }
#  }
#}
#dtplot1 =dtplot
# zt=zt1[pname%in%c("BSTT","BSTB","BSTR","TBS"),]
# zt$gev<- paste("GxE0.",zt$gev,sep = "")
# dtplot1=dtplot[dttype==2,]
# dtplot1$gev<- paste("GxE0.",dtplot1$gev,sep = "")
# for (p in c("BSTT","BSTB","BSTR","TBS")) {
#   for(g in c("20%","50%","80%","pblup")){
#     for(ge in c("GxE0.2","GxE0.5","GxE0.8")){
#       geneticgain= dtplot1[pname==p&genoratio==g&gev==ge,gv]-dtplot1[pname==p&genoratio==g&gev==ge,gv][1]
#       dtplot1[pname==p&genoratio==g&gev==ge,gg0:=geneticgain]
#       newgg= (dtplot1[pname==p&genoratio==g&gev==ge,gg0] - zt[pname==p&genoratio==g&gev==ge,gg0])/zt[pname==p&genoratio==g&gev==ge,gg0se]
#       dtplot1[pname==p&genoratio==g&gev==ge,gg0:=newgg]
#     }
#   }
# }
# dtplot1$pname = paste(dtplot1$pname,dtplot1$gev,sep="")
# picture<- ggplot(data = dtplot1,aes(generation,y=gg0,group=pname,color=gev))+
#   #geom_point()+
#   geom_line()+
#   xlab("Generation")+
#   #scale_shape_manual(values = c(3,15,0,NA))+
#   ylab("degree of deviation (standard variance)")+
#   theme_bw()+
#   theme(panel.grid.major = element_line(colour = NA),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank(),
#         text = element_text(family = "STXihei"),
#         legend.position = "right",
#         legend.background = element_rect(colour = "black"))+labs(color="Initial GEI")+
#   #geom_ribbon(aes(ymin = gxe -gxese, ymax = gxe+gxese), alpha = 0.3)+
#   scale_x_continuous(limits = c(0,30),breaks = seq(0,30,10))+facet_wrap(~ gev+genoratio)+ggtitle("The deviation from the mean genetic gain for all schemes\n when the GxE were the lowest one of 10 repeats") + theme(plot.title = element_text(hjust = 0.5))
# picture

dtplot1=dtplot[dttype==2,]
dtplot1[generation==0&gev==2,gxe:=0.2]
dtplot1[generation==0&gev==5,gxe:=0.5]
dtplot1[generation==0&gev==8,gxe:=0.8]
dtplot1[gev==2,gev:=0.2]
dtplot1[gev==5,gev:=0.5]
dtplot1[gev==8,gev:=0.8]
dtplot1$gev<- as.character(dtplot1$gev)
#dtplot1$gev<- paste("GxE0.",dtplot1$gev,sep = "")
dtplot1$pname = paste(dtplot1$pname,dtplot1$gev,sep="")
#dtplot1$pname = paste(dtplot1$pname,dtplot1$dttype,sep = "")
#dtplot1$dttype= as.character(dtplot1$dttype)
dtplot1[genoratio=="20%",genoratio:=20]
dtplot1[genoratio=="50%",genoratio:=50]
dtplot1[genoratio=="80%",genoratio:=80]

dtplot1[genoratio=="pblup",method:="Pedigree-based selection"]
dtplot1[genoratio!="pblup",method:="Genomic selection"]

dtplot1[genoratio=="pblup",genoratio:=0]

picture<- ggplot(data = dtplot1,aes(generation,y=gxe,group=pname,color=gev))+
  #geom_point()+
  geom_line()+
  xlab("Generation")+
  #scale_shape_manual(values = c(3,15,0,NA))+
  ylab("GEI intensity")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text = element_text(family = "STXihei"),
        legend.position = "right",
        legend.background = element_rect(colour = "black"))+#labs(color="Initial GEI")+
  #geom_ribbon(aes(ymin = gxe -gxese, ymax = gxe+gxese), alpha = 0.3)+
theme_zg()+labs(color = "GEI")+  #theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())+
  scale_x_continuous(limits = c(0,30),breaks = seq(0,30,10))+facet_wrap(~ method+genoratio)+scale_color_aaas()#+ggtitle("The GEI of biosecurity schemes during selective breeding of 30 generations\n (the lowest one of 10 repeats)") + theme(plot.title = element_text(hjust = 0.5))
picture

ggsave("gxe30th.png", picture, width = 10, height = 5, dpi = 300)
ggsave("gxeall.png", picture, width = 10, height = 5, dpi = 300)
ggsave("Figrue6.pdf", picture, width = 10, height = 5, dpi = 300)

# go="80%"
# ge=2
# dtplot[pname=="BSTB"&genoratio==go&gev==ge&dttype==2&generation==(which.min(dtplot[pname=="BSTB"&genoratio==go&gev==ge&dttype==2,gxe])-1),.(pname,generation,round(gxe,2))]
# dtplot[pname=="BSTT"&genoratio==go&gev==ge&dttype==2&generation==(which.min(dtplot[pname=="BSTT"&genoratio==go&gev==ge&dttype==2,gxe])-1),.(pname,generation,round(gxe,2))]
# dtplot[pname=="BSTR"&genoratio==go&gev==ge&dttype==2&generation==(which.min(dtplot[pname=="BSTR"&genoratio==go&gev==ge&dttype==2,gxe])-1),.(pname,generation,round(gxe,2))]
# 
# dtplot[pname=="TBS"&gev==ge&dttype==2&generation==(which.min(dtplot[pname=="TBS"&gev==ge&dttype==2,gxe])-1),.(pname,generation,round(gxe,2))]

# zt= zt1[pname%in%c("BSTT","BSTB","BSTR","TBS"),]
# zt$gev = paste("0.",zt$gev,sep="")
# 
# # zt$sem = zt$pname
# zt$pname = paste(zt$pname,zt$gev,sep="")
# 
# zt[genoratio=="20%",genoratio:=20]
# zt[genoratio=="50%",genoratio:=50]
# zt[genoratio=="80%",genoratio:=80]
# zt[genoratio=="pblup",method:="Pedigree-based selection"]
# zt[genoratio!="pblup",method:="Genomic selection"]
# 
# zt[genoratio=="pblup",genoratio:=0]
# 
# picture<- ggplot(data = zt,aes(nGeneration,y=gxxe,group=pname,fill=gev))+
#   # geom_point()+
#   geom_line()+
#   xlab("Generation")+
#   #scale_shape_manual(values = c(3,15,0,NA))+
#   ylab("GEI intensity")+
#   theme_bw()+
#   theme(panel.grid.major = element_line(colour = NA),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank(),
#         text = element_text(family = "STXihei"),
#         legend.position = "right",
#         legend.background = element_rect(colour = "black"))+labs(fill="GEI")+
#   geom_ribbon(aes(ymin = gxxe -gxese, ymax = gxxe+gxese), alpha = 0.3)+
#   #geom_errorbar(aes(ymin=gxxe-gxese,
#    #                 ymax=gxxe+gxese),
#     #            width=0.05,alpha = 0.5)+
#   scale_x_continuous(limits = c(0,30),breaks = seq(0,30,10))+facet_wrap(~ method+genoratio)+#ggtitle("The GEI of biosecurity schemes during selective breeding of 30 generations\n (average)") 
#   theme_zg()+scale_fill_aaas()
# picture
# ggsave("average.png", picture, width = 10, height = 5, dpi = 300)
# ggsave("Figrue5.pdf", picture, width = 10, height = 5, dpi = 300)
#dttype 1 = 30th
#dttype 2 = all
dtplot[dttype==1&generation==30&gev==2&genoratio=="20%",.(pname,gxe)]
zt1[nGeneration==30&gev==2&genoratio=="20%"&pname%in%c("BSTT","BSTB","BSTR"),.(pname,gxese)]

#se all

#gxe decreased
# zt<- zt1[pname%in%c("BSTT","BSTB","BSTR","TBS"),]
# zt$gxxe
# zt<- zt1[pname%in%c("TBS"),]
# zt$gev = as.character(zt$gev)
# zt[,pname:=paste(zt$pname,zt$gev,zt$genoratio,sep = "_")]
# picture<- ggplot(data = zt,aes(nGeneration,y=gxxe,group=pname,color=genoratio))+
#   geom_point()+
#   geom_line()+
#   xlab("Generation")+
#   #scale_shape_manual(values = c(3,15,0,NA))+
#   ylab("GxE interaction")+
#   theme_bw()+
#   theme(panel.grid.major = element_line(colour = NA),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank(),
#         text = element_text(family = "STXihei"),
#         legend.position = "right",
#         legend.background = element_rect(colour = "black"))+
#   scale_x_continuous(limits = c(0,30),breaks = seq(0,30,5))
# picture
# picture+ylim(0.1,0.85)
# lm(log(gxxe) ~ nGeneration ,data = zt[genoratio=="pblup"&gev=="2",])
zt = zt1[pname%in%c("BSTT","BSTB","BSTR","TBS"),]
zt$gev = paste("GEI 0.",zt$gev,sep="")
zt[pname == "BSTT",pname := "TOP"]
zt[pname == "BSTB",pname := "T&B"]
zt[pname == "BSTR",pname := "RAN"]
zt[pname == "TBS",pname := "PED"]
zt[genoratio=="20%",genoratio:=20]
zt[genoratio=="50%",genoratio:=50]
zt[genoratio=="80%",genoratio:=80]
zt[genoratio=="pblup",method:="Pedigree-based selection"]
zt[genoratio!="pblup",method:="Genomic selection"]
zt[genoratio=="pblup",genoratio:=0]

picture_gxe<- ggplot(data = zt,aes(nGeneration,y=gxxe,group=pname,color=pname))+
  #geom_point()+
  geom_line()+
  xlab("Generation")+
  #scale_shape_manual(values = c(3,15,0,NA))+
  ylab("GEI intensity")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text = element_text(family = "STXihei"),
        legend.position = "right",
        legend.background = element_rect(colour = "black"))+labs(color="Selective\ngenotyping\nmethods")+
  # geom_ribbon(aes(ymin = gg0 -gg0se, ymax = gg0+gg0se), alpha = 0.3)+
  geom_errorbar(aes(ymin=gxxe-gxese,
                    ymax=gxxe+gxese),
                width=0.05,alpha = 0.5)+
  scale_x_continuous(limits = c(0,30),breaks = seq(0,30,10))+facet_wrap(~ gev+genoratio)+#ggtitle("The GEI of biosecurity schemes during selective breeding of 30 generations\n (average)") 
  theme_zg()+scale_color_aaas()+theme(legend.title=element_blank())
picture_gxe
ggsave("FirstRvision_Figure4.pdf", picture_gxe , width = 10, height = 5, dpi = 300)


picture_gg0_bsg<- ggplot(data = zt,aes(nGeneration,y=gg0,group=pname,color=pname))+
  #geom_point()+
  geom_line()+
  xlab("Generation")+
  #scale_shape_manual(values = c(3,15,0,NA))+
  ylab("Genetic gain")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text = element_text(family = "STXihei"),
        legend.position = "right",
        legend.background = element_rect(colour = "black"))+labs(color="Selective\ngenotyping\nmethods")+
  # geom_ribbon(aes(ymin = gg0 -gg0se, ymax = gg0+gg0se), alpha = 0.3)+
  geom_errorbar(aes(ymin=gg0-gg0se,
                  ymax=gg0+gg0se),
             width=0.05,alpha = 0.5)+
  scale_x_continuous(limits = c(0,30),breaks = seq(0,30,10))+facet_wrap(~ gev+genoratio)+#ggtitle("The GEI of biosecurity schemes during selective breeding of 30 generations\n (average)") 
  theme_zg()+scale_color_aaas()+theme(legend.title=element_blank())
picture_gg0_bsg
ggsave("AD1_bsg.pdf", picture_gg0_bsg , width = 10, height = 5, dpi = 300)

picture_covig<- ggplot(data = zt,aes(nGeneration,y=covig,group=pname,color=pname))+
  #geom_point()+
  geom_line()+
  xlab("Generation")+
  #scale_shape_manual(values = c(3,15,0,NA))+
  ylab("Covariance")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text = element_text(family = "STXihei"),
        legend.position = "right",
        legend.background = element_rect(colour = "black"))+labs(color=NULL)+
  # geom_ribbon(aes(ymin = covig -covigse, ymax = covig+covigse), alpha = 0.3)+
  geom_errorbar(aes(ymin=covig -covigse,
                    ymax=covig +covigse),
                width=0.05,alpha = 0.5)+
  scale_x_continuous(limits = c(0,30),breaks = seq(0,30,10))+facet_wrap(~ gev+genoratio)+#ggtitle("The GEI of biosecurity schemes during selective breeding of 30 generations\n (average)") 
  theme_zg()+scale_color_aaas()
picture_covig
ggsave("AD3_cov.pdf", picture_covig , width = 10, height = 5, dpi = 300)

picture_i2<- ggplot(data = zt,aes(nGeneration,y=i2,group=pname,color=pname))+
  #geom_point()+
  geom_line()+
  xlab("Generation")+
  #scale_shape_manual(values = c(3,15,0,NA))+
  ylab("Interaction variance")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text = element_text(family = "STXihei"),
        legend.position = "right",
        legend.background = element_rect(colour = "black"))+labs(color=NULL)+
  # geom_ribbon(aes(ymin = covig -covigse, ymax = covig+covigse), alpha = 0.3)+
  geom_errorbar(aes(ymin=i2 -i2se,
                    ymax=i2 +i2se),
                width=0.05,alpha = 0.5)+
  scale_x_continuous(limits = c(0,30),breaks = seq(0,30,10))+facet_wrap(~ gev+genoratio)+#ggtitle("The GEI of biosecurity schemes during selective breeding of 30 generations\n (average)") 
  theme_zg()+scale_color_aaas()
picture_i2+ylim(0,5)
ggsave("AD3_I2.pdf", picture_i2+ylim(0,5) , width = 10, height = 5, dpi = 300)

picture_g2<- ggplot(data = zt,aes(nGeneration,y=g2,group=pname,color=pname))+
  #geom_point()+
  geom_line()+
  xlab("Generation")+
  #scale_shape_manual(values = c(3,15,0,NA))+
  ylab("Variance of the mean genetic effect")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text = element_text(family = "STXihei"),
        legend.position = "right",
        legend.background = element_rect(colour = "black"))+labs(color=NULL)+
  # geom_ribbon(aes(ymin = covig -covigse, ymax = covig+covigse), alpha = 0.3)+
  geom_errorbar(aes(ymin=g2 -g2se,
                    ymax=g2 +g2se),
                width=0.05,alpha = 0.5)+
  scale_x_continuous(limits = c(0,30),breaks = seq(0,30,10))+facet_wrap(~ gev+genoratio)+#ggtitle("The GEI of biosecurity schemes during selective breeding of 30 generations\n (average)") 
  theme_zg()+scale_color_aaas()
picture_g2+ylim(0,5)
ggsave("AD3_G2.pdf", picture_g2+ylim(0,5) , width = 10, height = 5, dpi = 300)

picture_idg<- ggplot(data = zt,aes(nGeneration,y=idg,group=pname,color=pname))+
  #geom_point()+
  geom_line()+
  xlab("Generation")+
  #scale_shape_manual(values = c(3,15,0,NA))+
  ylab("Ratio")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text = element_text(family = "STXihei"),
        legend.position = "right",
        legend.background = element_rect(colour = "black"))+labs(color=NULL)+
  # geom_ribbon(aes(ymin = covig -covigse, ymax = covig+covigse), alpha = 0.3)+
  geom_errorbar(aes(ymin=idg - idgse,
                    ymax=idg + idgse),
                width=0.05,alpha = 0.5)+
  scale_x_continuous(limits = c(0,30),breaks = seq(0,30,10))+facet_wrap(~ gev+genoratio)+#ggtitle("The GEI of biosecurity schemes during selective breeding of 30 generations\n (average)") 
  theme_zg()+scale_color_aaas()
picture_idg
ggsave("AD3_IG.pdf", picture_idg , width = 10, height = 5, dpi = 300)



zt = zt1[pname%in%c("TDS","DST"),]
zt$gev = paste("GEI 0.",zt$gev,sep="")
zt[genoratio=="pblup",method:="Pedigree-based selection"]
zt[genoratio!="pblup",method:="Genomic selection"]
picture_gg0_nbs<- ggplot(data = zt,aes(nGeneration,y=gg0,group=method,color=method))+
  #geom_point()+
  geom_line()+
  xlab("Generation")+
  #scale_shape_manual(values = c(3,15,0,NA))+
  ylab("Genetic gain")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text = element_text(family = "STXihei"),
        legend.position = "right",
        legend.background = element_rect(colour = "black"))+labs(color=NULL)+
  # geom_ribbon(aes(ymin = gg0 -gg0se, ymax = gg0+gg0se), alpha = 0.3)+
  geom_errorbar(aes(ymin=gg0-gg0se,
                  ymax=gg0+gg0se),
             width=0.05,alpha = 0.5)+
  scale_x_continuous(limits = c(0,30),breaks = seq(0,30,10))+facet_wrap(~ gev)+#ggtitle("The GEI of biosecurity schemes during selective breeding of 30 generations\n (average)") 
  theme_zg()+scale_color_aaas()
picture_gg0_nbs
ggsave("AD1_nbs.pdf", picture_gg0_nbs , width = 10, height = 5, dpi = 300)



ztforva = function(zt1,envtype){
  zt= zt1[pname%in%c("BSTT","BSTB","BSTR","TBS"),]
  zt$gev = paste("GEI 0.",zt$gev,sep="")
  zt[pname == "BSTT",pname := "TOP"]
  zt[pname == "BSTB",pname := "T&B"]
  zt[pname == "BSTR",pname := "RAN"]
  zt[pname == "TBS",pname := "PED"]
  zt[genoratio=="20%",genoratio:=20]
  zt[genoratio=="50%",genoratio:=50]
  zt[genoratio=="80%",genoratio:=80]
  zt[genoratio=="pblup",method:="Pedigree-based selection"]
  zt[genoratio!="pblup",method:="Genomic selection"]
  zt[genoratio=="pblup",genoratio:=0]
  if(envtype == "CE"){
    zt[,va1:=zt$va2]
    zt[,va1se:=zt$va2se]
    zt$envtype = "CE"
  }else{
    zt$envtype = "NE"
  }
  # zt$pname = paste(zt$pname,zt$envtype,sep = "")
  
  return(zt)
}

zt_CE = ztforva(zt1 = zt1, envtype = "CE")
zt_NE = ztforva(zt1 = zt1, envtype = "NE")
zt = rbind(zt_CE,zt_NE)

picture_va = ggplot(data = zt,aes(nGeneration,y=va1,group=envtype,color = pname,shape = envtype))+
  geom_point()+
  geom_line()+
  xlab("Generation")+
  #scale_shape_manual(values = c(3,15,0,NA))+
  ylab("Genetic variance")+
  theme_bw()+
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        text = element_text(family = "STXihei"),
        legend.position = "right",
        legend.background = element_rect(colour = "black"))+labs(color=NULL,shape = NULL)+
  # geom_ribbon(aes(ymin = covig -covigse, ymax = covig+covigse), alpha = 0.3)+
  geom_errorbar(aes(ymin=va1 -va1se,
                    ymax=va1 +va1se),
                width=0.05,alpha = 0.5)+
  scale_x_continuous(limits = c(0,30),breaks = seq(0,30,10))+facet_wrap(~ gev+genoratio)+#ggtitle("The GEI of biosecurity schemes during selective breeding of 30 generations\n (average)") 
  theme_zg()+scale_color_aaas()
picture_va
ggsave("AD2_va.pdf", picture_va , width = 10, height = 5, dpi = 300)



# zt= zt1[pname%in%c("BSTT","BSTB","BSTR","TBS"),]
# zt$gev = paste("0.",zt$gev,sep="")
# zt$pname = paste(zt$pname,zt$gev,sep="")
# 
# zt[genoratio=="20%",genoratio:=20]
# zt[genoratio=="50%",genoratio:=50]
# zt[genoratio=="80%",genoratio:=80]
# zt[genoratio=="pblup",method:="Pedigree-based selection"]
# zt[genoratio!="pblup",method:="Genomic selection"]
# 
# zt[genoratio=="pblup",genoratio:=0]


# picture_va1<- ggplot(data = zt,aes(nGeneration,y=va1,group=pname,color=gev))+
#   #geom_point()+
#   geom_line()+
#   xlab("Generation")+
#   #scale_shape_manual(values = c(3,15,0,NA))+
#   ylab("genetic variance in NE")+
#   theme_bw()+
#   theme(panel.grid.major = element_line(colour = NA),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank(),
#         text = element_text(family = "STXihei"),
#         legend.position = "right",
#         legend.background = element_rect(colour = "black"))+labs(color="GEI")+
#   # geom_ribbon(aes(ymin = va1 -va1se, ymax = va1+va1se), alpha = 0.3)+
#   geom_errorbar(aes(ymin= va1-va1se,
#                   ymax= va1+va1se),
#              width=0.05,alpha = 0.5)+
#   scale_x_continuous(limits = c(0,30),breaks = seq(0,30,10))+facet_wrap(~ method+genoratio)+#ggtitle("The GEI of biosecurity schemes during selective breeding of 30 generations\n (average)") 
#   theme_zg()+scale_color_aaas()+ylim(1,6)
# picture_va1
# 
# picture_va2<- ggplot(data = zt,aes(nGeneration,y=va2,group=pname,color=gev))+
#   #geom_point()+
#   geom_line()+
#   xlab("Generation")+
#   #scale_shape_manual(values = c(3,15,0,NA))+
#   ylab("genetic variance in CE")+
#   theme_bw()+
#   theme(panel.grid.major = element_line(colour = NA),
#         panel.background = element_rect(fill = "transparent",colour = NA),
#         plot.background = element_rect(fill = "transparent",colour = NA),
#         panel.grid.minor = element_blank(),
#         text = element_text(family = "STXihei"),
#         legend.position = "right",
#         legend.background = element_rect(colour = "black"))+labs(color="GEI")+
#   # geom_ribbon(aes(ymin = va2 -va2se, ymax = va2+va2se), alpha = 0.3)+
#   geom_errorbar(aes(ymin=va2-va2se,
#                   ymax=va2+va2se),
#              width=0.05,alpha = 0.5)+
#   scale_x_continuous(limits = c(0,30),breaks = seq(0,30,10))+facet_wrap(~ method+genoratio)+#ggtitle("The GEI of biosecurity schemes during selective breeding of 30 generations\n (average)") 
#   theme_zg()+scale_color_aaas()+ylim(1,6)
# picture_va2

# figure_va <- ggarrange(picture_va1 +rremove("ylab")+ylim(1,6) , picture_va2 + rremove("ylab")+ylim(1,6) , # remove axis labels from plots
#                        labels = c("a","b"),
#                        ncol = 2, nrow =1,
#                        widths = c(1,1),
#                        common.legend = TRUE, legend = "right",
#                        align = "v",
#                        font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
# figure_va
# ggsave("va12.png", figure_ig , width = 10, height = 5, dpi = 300)
# 
# 
# figure_ig <- ggarrange(picture_g2 +rremove("ylab")+ylim(0,5) , picture_i2 + rremove("ylab")+ylim(0,5) , # remove axis labels from plots
#                    labels = c("a","b"),
#                   ncol = 2, nrow =1,
#                    widths = c(1,1),
#                   common.legend = TRUE, legend = "right",
#                    align = "v",
#                   font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
# figure_ig
# ggsave("igvar.png", figure_ig , width = 10, height = 5, dpi = 300)