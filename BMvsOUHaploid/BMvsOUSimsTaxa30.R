rm(list=ls())
source("~/Dropbox/YiWeiHsu/Rcode/main.code/BMsnp/BMsnp.R")
source("~/Dropbox/YiWeiHsu/Rcode/main.code/OUsnp/OUsnp.R")

setwd("~/Dropbox/YiWeiHsu/Rcode/simulation/BMvsOUHaploid/Taxa30")

k.array <- 3:5#3:6
num.tree<-10
num.repdata<-50
seqlen<-1000
# ouparams<-c(5,5,40,40,90,80,100)
ouparams<-c(0.25,0.25,1,1,0,0,0)
names(ouparams)<-c("alp1","alp2","sig1","sig2","th0","th1","th2")
taxa.size.array<-c(30) #c(10,30,50)
sims<-500

obs.LSS.bm.array <- array(0,c(num.tree,num.repdata,length(k.array)))
obs.LSS.ou.array<- array(0,c(num.tree,num.repdata,length(k.array)))
sim.LSS.bm.array<-array(0,c(num.tree,num.repdata,length(k.array),sims))
sim.LSS.ou.array<- array(0,c(num.tree,num.repdata,length(k.array),sims))

for (taxaIndex in 1:length(taxa.size.array)){
  #taxaIndex <-1
  taxa.size<-taxa.size.array[taxaIndex]
  for (treeIndex in 1:num.tree){
    #treeIndex<-1
    #tree.trait<-sim.snp.tree.yg(taxa.size=taxa.size,seqlen=seqlen,ouparams=ouparams)
    tree<-sim.snp.tree(taxa.size=taxa.size,seqlen=seqlen)
    #plot(tree)
    tree<-chronos(tree)
    plot(tree)
    tree$tip.label<-paste("t",tree$tip.label,sep="")
    vcv.tree<-vcv(tree)
    diag(vcv.tree)<-max(vcv.tree)
    vcv.tree<-vcv.tree/max(vcv.tree)
    try(tree<-vcv2phylo(vcv.tree))
    plot(tree)
    
    #count_fail<-0
    for(repdataIndex in 1: num.repdata){
      #repdataIndex<-1
      trait.X<-sim.snp.yg(tree=tree,ouparams=ouparams)#trait only
      trait<-matrix(trait.X$X,ncol=1)
      rownames(trait)<-trait.X$Genus_species
      
      
      for(kIndex in 1:length(k.array)){
        #kIndex<-1
        DV.data<-NULL
        try(DV.data<-DV(k=k.array[kIndex],tree=tree))
        if(!is.null(DV.data)){
          V<-DV.data$V
          D<-DV.data$D
          
          obs.lss.bm<-NULL
          try(obs.lss.bm <- LSS.bm(k=k.array[kIndex], trait=trait, V=V,D=D))
          if(!is.null(obs.lss.bm)){
            obs.LSS.bm.array[treeIndex,repdataIndex,kIndex] <- obs.lss.bm$lss
            for(simIndex in 1:sims){
              print(paste("taxa=", taxa.size.array[taxaIndex], ", tree=", treeIndex, ", repdata=", repdataIndex, ", k=", k.array[kIndex] , ", sim=", simIndex , sep=""))
              #simIndex<-1
              trait<-sample(trait)
              sim.lss.bm<-NULL
              try(sim.lss.bm<-LSS.bm(k=k.array[kIndex],trait=trait,V=V,D=D))
              if(!is.null(sim.lss.bm)){
                try(sim.LSS.bm.array[treeIndex,repdataIndex,kIndex,simIndex]<-sim.lss.bm$lss)
              }else{sim.LSS.bm.array[treeIndex,repdataIndex,kIndex,simIndex]<--Inf}
            }#end of bm sims
          }else{
            obs.LSS.bm.array[treeIndex,repdataIndex,kIndex]<-  -Inf
            sim.LSS.bm.array[treeIndex,repdataIndex,kIndex,]<- -Inf
          }#end of bm is.null
          
          obs.lss.ou<-NULL
          try(obs.lss.ou <- LSS.ou(k=k.array[kIndex], trait=trait, V=V, D=D))
          if(!is.null(obs.lss.ou)){
            obs.LSS.ou.array[treeIndex, repdataIndex,kIndex] <- obs.lss.ou$lss
            for(simIndex in 1:sims){
              #simIndex<-1
              trait<-sample(trait)
              sim.lss.ou<-NULL
              try(sim.lss.ou<-LSS.ou(k=k.array[kIndex],trait=trait,V=V,D=D))
              if(!is.null(sim.lss.ou)){
                try(sim.LSS.ou.array[treeIndex, repdataIndex,kIndex,simIndex]<-sim.lss.ou$lss)
              }else{sim.LSS.ou.array[treeIndex, repdataIndex,kIndex,simIndex]<--Inf}
            }#end of ou sims
          }else{
            obs.LSS.ou.array[treeIndex, repdataIndex,kIndex]<-  -Inf
            sim.LSS.ou.array[treeIndex, repdataIndex,kIndex,]<- -Inf
          }#end of ou is.null
        }else{
          #print(paste( "dataIndex=",repdataIndex,", failed k=",k.array[kIndex],sep=""))
          #count_fail<-count_fail+1
          obs.LSS.bm.array[treeIndex,repdataIndex,kIndex]<-  -Inf
          sim.LSS.bm.array[treeIndex,repdataIndex,kIndex,]<- -Inf
          obs.LSS.ou.array[treeIndex,repdataIndex,kIndex]<-  -Inf
          sim.LSS.ou.array[treeIndex,repdataIndex,kIndex,]<- -Inf
        }#end of DV
      }#end of kIndex
      
      
      max.LSS.bm<-max(obs.LSS.bm.array[treeIndex,repdataIndex,])
      max.LSS.ou<-max(obs.LSS.ou.array[treeIndex,repdataIndex,])
      
      best.k.bm<-unique(k.array[which(obs.LSS.bm.array[treeIndex,repdataIndex,]==max.LSS.bm)])[1]
      best.k.ou<-unique(k.array[which(obs.LSS.ou.array[treeIndex,repdataIndex,]==max.LSS.ou)])[1]
      
      max.sim.LSS.bm.array<-apply(sim.LSS.bm.array[treeIndex,repdataIndex,,],2,max)
      max.sim.LSS.ou.array<-apply(sim.LSS.ou.array[treeIndex,repdataIndex,,],2,max)
      
      p.value.bm<-sum(max.sim.LSS.bm.array>c(max.LSS.bm))/sims
      p.value.ou<-sum(max.sim.LSS.ou.array>c(max.LSS.ou))/sims
      print(paste("bm p-value=", p.value.bm ,  ", ou p-value=" , p.value.ou, sep=""))
      
      
      bm.DV.data<-NULL
      try(bm.DV.data<-DV(k=best.k.bm,tree=tree))
      obs.lss.best.k.bm<-NULL
      try(obs.lss.best.k.bm<-LSS.bm(k=best.k.bm, trait=trait, V=bm.DV.data$V,D=bm.DV.data$D))
      
      ou.DV.data<-NULL
      try(ou.DV.data<-DV(k=best.k.ou,tree=tree))
      obs.lss.best.k.ou<-NULL
      try(obs.lss.best.k.ou<-LSS.ou(k=best.k.ou, trait=trait, V=ou.DV.data$V,D=ou.DV.data$D))
      
      output<-list(tree=tree, trait=trait, best.k.bm=best.k.bm, best.k.ou=best.k.ou, bm.DV.data=bm.DV.data, ou.DV.data=ou.DV.data, obs.lss.best.k.bm=obs.lss.best.k.bm,  obs.LSS.bm.array=obs.LSS.bm.array, sim.LSS.bm.array=sim.LSS.bm.array, obs.lss.best.k.ou=obs.lss.best.k.ou, obs.LSS.ou.array=obs.LSS.ou.array, sim.LSS.ou.array=sim.LSS.ou.array, p.value.bm=p.value.bm, p.value.ou=p.value.ou)
      assign(paste("taxa" , taxa.size , "tree" , treeIndex, "data",repdataIndex , sep=""),output)
      setwd("/Users/yiweihsu/Documents/BMvsOUSimsHaploid/")
      save.image(file="BMvsOUSimsTaxa30.RData")
      setwd("~/Dropbox/YiWeiHsu/Rcode/simulation/BMvsOUHaploid/Taxa30")
      
    }#end of repdataIndex
  }#end of treeIndex
}#end of taxaIndex



