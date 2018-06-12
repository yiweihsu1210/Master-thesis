rm(list = ls())
setwd("/Users/yiweihsu/Documents/OUsnpPowerv2")
system("open .")
system("pwd")
system("ls")

library(xtable)

RDfile<-c("OU16.RData","OU32.RData","OU64.RData","OU128.RData")
#oupow<-array(NA,c(6,4))
#bmpow<-array(NA,c(6,4))

bm.power.alpha1.sig.array<-NULL
bm.power.alpha2.sig.array<-NULL
ou.power.alpha1.sig.array<-NULL
ou.power.alpha2.sig.array<-NULL


alp1 <- alp2 <-c(5, 10) 
sig1 <-sig2 <- c(1, 5, 10)
taxa.size.array <- c(16,32,64,128) 
num.tree<-5
num.repdata<-100


p.value.array.bm<-array(0,c(length(RDfile),length(alp1),length(sig1),length(taxa.size.array),num.tree,num.repdata ))
p.value.array.ou<-array(0,c(length(RDfile),length(alp1),length(sig1),length(taxa.size.array),num.tree,num.repdata ))



for(fileIndex in 1:length(RDfile)){
#  fileIndex<-1
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
  
  
  print("BM")
  #print(apply(power.bm.array[1,1,1,1],1,mean))
  
  power.bm.array
  
  bm.power.alpha1.sig.array<-rbind(bm.power.alpha1.sig.array, apply(power.bm.array[1,,1,],1,mean))
  bm.power.alpha2.sig.array<-rbind(bm.power.alpha2.sig.array, apply(power.bm.array[2,,1,],1,mean))
  
  
  print(apply(power.bm.array[2,,1,],1,mean))
  
  print("OU")
  print(apply(power.ou.array[1,,1,],1,mean))
  print(apply(power.ou.array[2,,1,],1,mean))
  ou.power.alpha1.sig.array<-rbind(ou.power.alpha1.sig.array, apply(power.ou.array[1,,1,],1,mean))
  ou.power.alpha2.sig.array<-rbind(ou.power.alpha2.sig.array, apply(power.ou.array[2,,1,],1,mean))
  
  
  # 
  #  bmpow[,fileIndex]<-c(power.bm.array)
  #  oupow[,fileIndex]<-c(power.ou.array)
}
#bmpow
#oupow

#plot(power.bm.array)
#plot(power.ou.array)

print("OUdata")
colnames(bm.power.alpha1.sig.array)<-colnames(bm.power.alpha2.sig.array)<-sig1
colnames(ou.power.alpha1.sig.array)<-colnames(ou.power.alpha2.sig.array)<-sig1
rownames(bm.power.alpha1.sig.array)<-rownames(bm.power.alpha2.sig.array)<-c("16","32","64","128")
rownames(ou.power.alpha1.sig.array)<-rownames(ou.power.alpha2.sig.array)<-c("16","32","64","128")
print("BM")
print(alp1[1])
bm.power.alpha1.sig.array<- bm.power.alpha1.sig.array/num.repdata
print(alp1[2])
bm.power.alpha2.sig.array<-bm.power.alpha2.sig.array/num.repdata
print("OU")
print(alp1[1])
ou.power.alpha1.sig.array<-ou.power.alpha1.sig.array/num.repdata
print(alp1[2])
ou.power.alpha2.sig.array<-ou.power.alpha2.sig.array/num.repdata


xtable(bm.power.alpha1.sig.array)
xtable(bm.power.alpha2.sig.array)

xtable(ou.power.alpha1.sig.array)
xtable(ou.power.alpha2.sig.array)



#do plot
#alpIndex<-1
#sigIndex<-1
#taxaIndex<-1

#fileIndex<-1
#plot(sort(apply(p.value.array.bm[fileIndex,alpIndex,sigIndex,taxaIndex,,],2,mean,na.rm=T)))
#fileIndex<-2
#points(sort(apply(p.value.array.bm[fileIndex,alpIndex,sigIndex,taxaIndex,,],2,mean,na.rm=T)),col="red")
#fileIndex<-3
#points(sort(apply(p.value.array.bm[fileIndex,alpIndex,sigIndex,taxaIndex,,],2,mean,na.rm=T)),col="blue")
#fileIndex<-4
#points(sort(apply(p.value.array.bm[fileIndex,alpIndex,sigIndex,taxaIndex,,],2,mean,na.rm=T)),col="purple")

treeIndex<-1

par(mfrow=c(1,2))
for(alpIndex in c(1,2)){
  #  taxa16<-apply(p.value.array.bm[1,alpIndex,sigIndex,taxaIndex,,],2,mean,na.rm=T)
  #  taxa32<-apply(p.value.array.bm[2,alpIndex,sigIndex,taxaIndex,,],2,mean,na.rm=T)
  #  taxa64<-apply(p.value.array.bm[3,alpIndex,sigIndex,taxaIndex,,],2,mean,na.rm=T)
  #  taxa128<-apply(p.value.array.bm[4,alpIndex,sigIndex,taxaIndex,,],2,mean,na.rm=T)
  #  boxplot(taxa16,taxa32,taxa64,taxa128,main=paste("BM sigma = ",sig1[sigIndex],sep=""),names=c("16","32","64","128"),xlab="taxa size",ylab="p-value",ylim=c(0,0.8))
  dim(p.value.array.bm[1,alpIndex, ,taxaIndex,treeIndex,])
  
  taxa16<-apply(p.value.array.bm[1,alpIndex, ,taxaIndex,treeIndex,],2,mean,na.rm=T)
  taxa32<-apply(p.value.array.bm[2,alpIndex, ,taxaIndex,treeIndex,],2,mean,na.rm=T)
  taxa64<-apply(p.value.array.bm[3,alpIndex, ,taxaIndex,treeIndex,],2,mean,na.rm=T)
  taxa128<-apply(p.value.array.ou[4,alpIndex, ,taxaIndex,treeIndex,],2,mean,na.rm=T)
#  boxplot(taxa16,taxa32,taxa64,taxa128,main=paste("BM Data vs OU model: sigma = ",sig1[sigIndex],sep=""),names=c("16","32","64","128"),xlab="taxa size",ylab="p-value",ylim=c(0,1))
  boxplot(taxa16,taxa32,taxa64,taxa128,main=paste("OU Data vs BM model: alpha = ", alp1[alpIndex],sep=""),names=c("16","32","64","128"),xlab="taxa size",ylab="p-value",ylim=c(0,1))
  }
#p-value does not change via sigma, #p-value increases with sample size for BM it makes sense as data
#come from the BM model.

#lets look at ou. BM data so it is expected that BM has higher p-value(lower rejection rate given significant level)





