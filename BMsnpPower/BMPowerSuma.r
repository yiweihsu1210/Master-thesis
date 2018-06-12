rm(list = ls())
setwd("/Users/yiweihsu/Documents/BMsnpPowerv2")
system("pwd")

library(xtable)

RDfile<-c("BM16.RData","BM32.RData","BM64.RData","BM128.RData")
#oupow<-array(NA,c(4,4))
#bmpow<-array(NA,c(4,4))
bm.power.sig.array<-NULL
ou.power.sig.array<-NULL

alp1 <- alp2 <- c(1e-6)
sig1 <-sig2 <- c(1,5,10)
taxa.size.array <- c(16,32,64,128) 
num.tree<-5
num.repdata<-100

p.value.array.bm<-array(0,c(length(RDfile),length(alp1),length(sig1),length(taxa.size.array),num.tree,num.repdata ))
p.value.array.ou<-array(0,c(length(RDfile),length(alp1),length(sig1),length(taxa.size.array),num.tree,num.repdata ))


for(fileIndex in 1:length(RDfile)){
  #fileIndex<-1
  print(RDfile[fileIndex])
  load(RDfile[fileIndex])
  
  for(alpIndex in 1:length(alp1)){
    #alpIndex<-1
    alp.value<-alp1[alpIndex]
    for(sigIndex in 1:length(sig1)){
      sig.value<-sig1[sigIndex]
      for (taxaIndex in 1:length(taxa.size.array)){
        #taxaIndex <-1
        taxa.size<-taxa.size.array[taxaIndex]
        for (treeIndex in 1:num.tree){
          #num.tree<-1
          for(repdataIndex in 1: num.repdata){
            #num.repdata<-4
            output<-try(get(paste("alp", alp.value, "sig", sig.value, "taxa" , taxa.size , "tree" , treeIndex, "data",repdataIndex ,sep="")))
            p.value.array.bm[fileIndex,alpIndex,sigIndex,taxaIndex,treeIndex,repdataIndex]<-output$p.value.bm
            p.value.array.ou[fileIndex,alpIndex,sigIndex,taxaIndex,treeIndex,repdataIndex]<-output$p.value.ou
          }
        }
      }
    }
  }
  
  
  
  
  
  
  print(RDfile[fileIndex])
  print("BM")
  print(apply(power.bm.array[1,,1,],1,mean))
  bm.power.sig.array<-rbind(bm.power.sig.array, apply(power.bm.array[1,,1,],1,mean))
  
  print("OU")
  print(apply(power.ou.array[1,,1,],1,mean))
  ou.power.sig.array<-rbind(ou.power.sig.array, apply(power.ou.array[1,,1,],1,mean))
  
  #bmpow[,fileIndex]<-c(power.bm.array)
  #oupow[,fileIndex]<-c(power.ou.array)
}
#bmpow
#oupow

#plot(power.bm.array)
#plot(power.ou.array)
print("BMdata")
colnames(bm.power.sig.array)<-colnames(ou.power.sig.array)<-sig1
rownames(bm.power.sig.array)<-rownames(ou.power.sig.array)<-c("16","32","64","128")

bm.power.sig.array<-bm.power.sig.array/num.repdata
ou.power.sig.array<-ou.power.sig.array/num.repdata

xtable(bm.power.sig.array,digits=2)
xtable(ou.power.sig.array,digits=2)

#Now plot p-value 
# coniditioned on tree
# length(RDfile) = 4
# length(alp1) = 1
# length(sig1) = 3
# length(taxa.size.array)= 4
# num.tree=5
# num.repdata=100
# dim(p.value.array.bm)
# 4   1   3   4   5 100
alpIndex<-1
sigIndex<-1
taxaIndex<-1

fileIndex<-1
plot(sort(apply(p.value.array.bm[fileIndex,alpIndex,sigIndex,taxaIndex,,],2,mean)))
fileIndex<-2
points(sort(apply(p.value.array.bm[fileIndex,alpIndex,sigIndex,taxaIndex,,],2,mean)),col="red")
fileIndex<-3
points(sort(apply(p.value.array.bm[fileIndex,alpIndex,sigIndex,taxaIndex,,],2,mean,na.rm=T)),col="blue")
fileIndex<-4
points(sort(apply(p.value.array.bm[fileIndex,alpIndex,sigIndex,taxaIndex,,],2,mean,na.rm=T)),col="purple")

par(mfrow=c(1,2))
for(sigIndex in c(1,2)){
#  taxa16<-apply(p.value.array.bm[1,alpIndex,sigIndex,taxaIndex,,],2,mean,na.rm=T)
#  taxa32<-apply(p.value.array.bm[2,alpIndex,sigIndex,taxaIndex,,],2,mean,na.rm=T)
#  taxa64<-apply(p.value.array.bm[3,alpIndex,sigIndex,taxaIndex,,],2,mean,na.rm=T)
#  taxa128<-apply(p.value.array.bm[4,alpIndex,sigIndex,taxaIndex,,],2,mean,na.rm=T)
#  boxplot(taxa16,taxa32,taxa64,taxa128,main=paste("BM sigma = ",sig1[sigIndex],sep=""),names=c("16","32","64","128"),xlab="taxa size",ylab="p-value",ylim=c(0,0.8))
  
  taxa16<-apply(p.value.array.ou[1,alpIndex,sigIndex,taxaIndex,,],2,mean,na.rm=T)
  taxa32<-apply(p.value.array.ou[2,alpIndex,sigIndex,taxaIndex,,],2,mean,na.rm=T)
  taxa64<-apply(p.value.array.ou[3,alpIndex,sigIndex,taxaIndex,,],2,mean,na.rm=T)
  taxa128<-apply(p.value.array.ou[4,alpIndex,sigIndex,taxaIndex,,],2,mean,na.rm=T)
  boxplot(taxa16,taxa32,taxa64,taxa128,main=paste("BM Data vs OU model: sigma = ",sig1[sigIndex],sep=""),names=c("16","32","64","128"),xlab="taxa size",ylab="p-value",ylim=c(0,1))
}
#p-value does not change via sigma, #p-value increases with sample size for BM it makes sense as data
#come from the BM model.

#lets look at ou. BM data so it is expected that BM has higher p-value(lower rejection rate given significant level)
