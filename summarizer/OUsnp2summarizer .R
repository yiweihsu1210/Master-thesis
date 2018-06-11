rm(list=ls())
setwd("/Users/yiweihsu/Documents/OUsnp2/")
load("/Users/yiweihsu/Documents/OUsnp2/OUdiploidSimsTaxa20.RData")
load("/Users/yiweihsu/Documents/OUsnp2/OUdiploidSimsTaxa30.RData")
load("/Users/yiweihsu/Documents/OUsnp2/OUdiploidSimsTaxa50.RData")

save.image("simoutest1.RData")

load("/Users/yiweihsu/Documents/OUsnp2/simoutest1.RData")

taxa.size.array<-c(20,30,50)
num.tree<-10
num.repdata<-50
rho.array<-seq(0,1,0.2)
k.array<- 3:5
sim.chains<-5000#50000


output.arrayk1<-array(0,c(length(taxa.size.array),num.tree,num.repdata,length(rho.array),sim.chains,8))
output.arrayk2<-array(0,c(length(taxa.size.array),num.tree,num.repdata,length(rho.array),sim.chains,9))
output.arrayk3<-array(0,c(length(taxa.size.array),num.tree,num.repdata,length(rho.array),sim.chains,10))
existed.count<-0
total.count<-0
for (taxaIndex in 1:length(taxa.size.array)) {
  #taxaIndex<-3
  #  load(paste("BMdiploidSimsTaxa",taxa.size.array[taxaIndex],".RData" , sep=""))
  #  taxa.size.array<-c(20,30,50)
  taxa.size<-taxa.size.array[taxaIndex]
  for (treeIndex in 1:num.tree) {
    #treeIndex<-1
    for (repdataIndex in 1:num.repdata) {
      #redataIndex<-10
      for(rhoIndex in 1:length(rho.array)){
       # rhoIndex<-2
        for (kIndex in 1:length(k.array)) {
          #kIndex<-1
          total.count<-total.count+1
          print(paste("taxa" , taxa.size , "tree" , treeIndex, "data",repdataIndex,"rho", rho.array[rhoIndex], "k", kIndex,sep=""))
          
          
          if(kIndex==1){
            object2get<-paste("taxa" , taxa.size , "tree" , treeIndex, "data",repdataIndex,"rho",  rho.array[rhoIndex], "k", kIndex,sep="")
            if(exists(object2get)){
              existed.count<-existed.count+1
              try(output <- get(object2get))
              if(!is.null(output)){
                output.arrayk1[taxaIndex,treeIndex,repdataIndex,rhoIndex,,]<-output$mhgibb.sample
                print(output$mhgibb.sample[1,])
              }
            }
          }
          
          if(kIndex==2){
            #gibb.sample<-array(0,c(sim.chains,8))
            #output<-list(gibb.sample=gibb.sample)
            object2get<-paste("taxa" , taxa.size , "tree" , treeIndex, "data",repdataIndex,"rho",  rho.array[rhoIndex], "k", kIndex,sep="")
            if(exists(object2get)){
              existed.count<-existed.count+1
              try(output <- get(object2get))
              if(!is.null(output)){
                output.arrayk2[taxaIndex,treeIndex,repdataIndex,rhoIndex,,]<-output$mhgibb.sample
                print(output$mhgibb.sample[1,])
              }
            }
          }
          if(kIndex==3){
            #gibb.sample<-array(0,c(sim.chains,9))
            #output<-list(gibb.sample=gibb.sample)
            object2get<-paste("taxa" , taxa.size , "tree" , treeIndex, "data",repdataIndex,"rho",  rho.array[rhoIndex], "k", kIndex,sep="")
            if(exists(object2get)){
              existed.count<-existed.count+1
              try(output <- get(object2get))
              if(!is.null(output)){
                output.arrayk3[taxaIndex,treeIndex,repdataIndex,rhoIndex,,]<-output$mhgibb.sample
                print(output$mhgibb.sample[1,])
              }
            }  
          }
          
          print(paste("finished",existed.count/total.count,sep=""))
          #dim(output.array[taxaIndex,treeIndex,repdataIndex,rhoIndex,kIndex,,])
          #dim(output$gibb.sample)
          #p.value.bm.array[repdataIndex]<-print(output$p.value.bm)
          #p.value.ou.array[repdataIndex]<-print(output$p.value.ou)
        }
      }
    }
  }
}
output.arrayk3



