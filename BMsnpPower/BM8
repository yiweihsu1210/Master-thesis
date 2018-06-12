rm(list=ls())
source("~/Dropbox/YiWeiHsu/Rcode/main.code/BMsnp/BMsnp.R")
source("~/Dropbox/YiWeiHsu/Rcode/main.code/OUsnp/OUsnp.R")

setwd("/Users/yiweihsu/Documents/BMsnpPowerv2")
system("pwd")

sign.level<-0.1
k.array <- c(3:5)#3:10
num.tree<- 5
num.repdata <- 100#200
seqlen<-1000
taxa.size.array <- c(16)  # c(50)  16,32,64,128
sims<- 100

# ouparams<-c(5,5,40,40,90,80,100)

alp1 <- alp2 <- c(1e-6)
sig1 <-sig2 <- c(1,5,10)#c(5, 20, 30, 40)
th0 <- 90
th1 <-80
th2 <-100

power.bm.array <- array(0,c(length(alp1),length(sig1),length(taxa.size.array),num.tree ))
power.ou.array <- array(0,c(length(alp1),length(sig1),length(taxa.size.array),num.tree ))

#dim(power.bm.array)



# library(TreeSim)
# tree <- sim.bd.taxa(n=taxa.size, numbsim=1,lambda=100,mu=0.5)[[1]]
# plot(tree)
# vcv.tree<-vcv(tree)
# diag(vcv.tree)<-max(vcv.tree)
# head(vcv.tree)

for(alpIndex in 1:length(alp1)){
  #alpIndex<-1
  alp.value<-alp1[alpIndex]
  for(sigIndex in 1:length(sig1)){
    #sigIndex<-1
    sig.value<-sig1[sigIndex]
    ouparams<-c(alp.value,alp.value,sig.value,sig.value,th0,th1,th2)
    names(ouparams)<-c("alp1","alp2","sig1","sig2","th0","th1","th2")
    print(ouparams)

    obs.LSS.bm.array <- array(0,c(num.tree,num.repdata,length(k.array)))
    obs.LSS.ou.array<- array(0,c(num.tree,num.repdata,length(k.array)))
    dim(obs.LSS.bm.array)
    sim.LSS.bm.array<-array(0,c(num.tree,num.repdata,length(k.array),sims))
    sim.LSS.ou.array<- array(0,c(num.tree,num.repdata,length(k.array),sims))

    for (taxaIndex in 1:length(taxa.size.array)){
      #taxaIndex <-1
      taxa.size<-taxa.size.array[taxaIndex]
      for (treeIndex in 1:num.tree){
        #treeIndex<-1
        tree.trait<-sim.snp.tree.yg(taxa.size=taxa.size,seqlen=seqlen,ouparams=ouparams)
        tree<-sim.snp.tree(taxa.size=taxa.size,seqlen=seqlen)
        plot(tree)
        tree<-chronos(tree)
        tree$tip.label<-paste("t",tree$tip.label,sep="")
        plot(tree)

        # tree<-stree(n=taxa.size, type = "balanced")
        # tree<-compute.brlen(tree)
        vcv.tree<-vcv(tree)
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
                  #print(paste("taxa=", taxa.size.array[taxaIndex], ", tree=", treeIndex, ", repdata=", repdataIndex, ", k=", k.array[kIndex] , ", sim=", simIndex , sep=""))
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

          
          
          #p.value.bm<-sum(max.sim.LSS.bm.array>c(max.LSS.bm))/sims
          #p.value.ou<-sum(max.sim.LSS.ou.array>c(max.LSS.ou))/sims
          
          p.value.bm<-NA
          p.value.ou<-NA
          
          if(max.LSS.bm!=-Inf){
          bmbdd <- max.LSS.bm + qt(p=sign.level,df=length(taxa.size-1),lower=FALSE)*sd(max.sim.LSS.bm.array)/sqrt(taxa.size)
          p.value.bm<-sum(max.sim.LSS.bm.array>c(bmbdd))/sims
          if(p.value.bm < sign.level){power.bm.array[alpIndex,sigIndex,taxaIndex,treeIndex] <- power.bm.array[alpIndex,sigIndex,taxaIndex,treeIndex] + 1}
          }
          
          if(max.LSS.ou!=-Inf){
          oubdd <- max.LSS.ou + qt(p=sign.level,df=length(taxa.size-1),lower=FALSE)*sd(max.sim.LSS.ou.array)/sqrt(taxa.size)
          p.value.ou<-sum(max.sim.LSS.ou.array>c(oubdd))/sims
          if(p.value.ou < sign.level){power.ou.array[alpIndex,sigIndex,taxaIndex,treeIndex] <- power.ou.array[alpIndex,sigIndex,taxaIndex,treeIndex] + 1}
          }
          print(paste("bm p-value=", p.value.bm ,  ", ou p-value=" , p.value.ou, sep=""))


          bm.DV.data<-NULL
          try(bm.DV.data<-DV(k=best.k.bm,tree=tree))
          obs.lss.best.k.bm<-NULL
          try(obs.lss.best.k.bm<-LSS.bm(k=best.k.bm, trait=trait, V=bm.DV.data$V,D=bm.DV.data$D))

          ou.DV.data<-NULL
          try(ou.DV.data<-DV(k=best.k.ou,tree=tree))
          obs.lss.best.k.ou<-NULL
          try(obs.lss.best.k.ou<-LSS.ou(k=best.k.ou, trait=trait, V=ou.DV.data$V,D=ou.DV.data$D))

          output<-list(p.value.bm=p.value.bm, p.value.ou=p.value.ou,power.bm.array=power.bm.array,power.ou.array=power.ou.array )
          #print("powerou BMdata")
          #print(power.ou.array)
          assign(paste("alp", alp.value, "sig", sig.value, "taxa" , taxa.size , "tree" , treeIndex, "data",repdataIndex ,sep=""),output)

          
          
          save.image(file="BM8.RData")


        }#end of repdataIndex
      }#end of treeIndex
    }#end of taxaIndex
    print("power bm")
    print(power.bm.array[alpIndex,sigIndex,taxaIndex,treeIndex]/num.repdata)
    print("power ou")
    print(power.ou.array[alpIndex,sigIndex,taxaIndex,treeIndex]/num.repdata)
  }#end sig.value
}#end.alp.value
