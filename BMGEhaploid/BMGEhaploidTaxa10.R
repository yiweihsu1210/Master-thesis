#move back to maincode folder after implementation
rm(list=ls())
source("~/Dropbox/YiWeiHsu/Rcode/main.code/BMsnp2_haploid/BMsnp2.r")
setwd("~/Dropbox/YiWeiHsu/Rcode/simulation/BMsnpGEhaploid/Taxa10/")

rho.array<-seq(0,1,0.2)#array(c(0.1),c(1)) #array(c(0.03162278),c(1))
#diploid.size.array <- c(10)
taxa.size.array <- c(10)# 2*diploid.size.array
k.array<- 3:4
r<-1
seqlen<-1000
ouparams<-c(0.5,0.5,5,5,90,80,100)
names(ouparams)<-c("alp1","alp2","sig1","sig2","th0","th1","th2")
num.tree<-5#10
num.repdata<-1#50
sim.chains<-50000
true.sigma.sq <- 5
true.b0=20
true.b1=40
true.nu.sq=2
true.beta<-c(true.b0,true.b1)

for (taxaIndex in 1:length(taxa.size.array)) {
  #taxaIndex <-1
  taxa.size<-taxa.size.array[taxaIndex]
  #generate a tree 
  #taxa.size<-taxa.size.array[1]
  for (treeIndex in 1:num.tree) {
    #treeIndex<-1
    tree<-sim.snp.tree(taxa.size=taxa.size,seqlen=seqlen)
    tree<-chronos(tree)
    #plot(tree)
    tree$tip.label<-paste("t",tree$tip.label,sep="")
    vcv.tree<-vcv(tree)
    diag(vcv.tree)<-max(vcv.tree)
    vcv.tree<-vcv.tree/max(vcv.tree)
    try(tree<-vcv2phylo(vcv.tree))
    #plot(tree)
    
    for(repdataIndex in 1:num.repdata){
      #repdataIndex<-1
      for(rhoIndex in 1:length(rho.array)){
        #rhoIndex<-1
        rho<-rho.array[rhoIndex]
        #generate trait : (1) genetic component (2) environmental component
        Zg.X <- sim.snp.yg(tree=tree,ouparams=ouparams)#trait only
        #Z<-d.mtx(n=diploid.size.array[1])
        #Zg<-Z%*%matrix(Zg.X$X,ncol=1)
        Zg<-matrix(Zg.X$X,ncol=1)
        
        #rownames(Zg) <- Zg.X$Genus_species
        Yg<-rho*Zg
        
        X<- matrix(runif(n=taxa.size,min=-5,max=5),ncol=1)
        
        
        #####
        #CHECK WE MAY NEED TO CHANGE 20, 4, 40 to other values for simulation.......
        #####
        Ze<-matrix(rmvnorm(n=1,mean=true.b0+true.b1*X,sigma=true.nu.sq*diag(1,c(taxa.size, taxa.size))),ncol=1)# may use differnt a,b,s for  a + bX, s*diag...
        Ye<-(1-rho)*Ze
        
        Y <- Yg+Ye
        
        #get clustering here
        for(kIndex in 1:length(k.array)){
          #kIndex<-1
          DV.data<-NULL
          try(DV.data<-DV(k=k.array[kIndex],tree=tree))
          if(!is.null(DV.data)){
            V<-DV.data$V
            D<-DV.data$D
            
            output<-NULL
            try(output<-Gibbs(Yg=Yg,Ye=Ye, k=k.array[kIndex],V=V, D=D,X=X, rho=rho,r=r, true.beta=true.beta, sigma.sq = true.sigma.sq, nu.sq = true.nu.sq , sim.chains=sim.chains))
            
            
            assign(paste("taxa" , taxa.size , "tree" , treeIndex, "data",repdataIndex,"rho", rho, "k", kIndex,sep=""),output)  
            #boxplot(output$gibb.sample[,2], main= "beta1")
            #boxplot(output$gibb.sample[,3], main="sigma.sq")
            #boxplot(output$gibb.sample[,4], main="nu.sq")
            
          }#end of DVisNULL
        }#end of kIndex
      }#end of rhoIndex
      setwd("/Users/yiweihsu/Documents/BMsnp2v2/")
      save.image(file="BMGEhaploidSimsTaxa10.RData")
      setwd("~/Dropbox/YiWeiHsu/Rcode/simulation/BMsnpGEhaploid/Taxa10/")
    }#end of repdataIndex
  }#end of treeIndex
}#end of taxaIndex

# 
# #return posterior from gibbs sampling
# par(mfrow=c(1,3))
# qnt<-quantile(output$gibb.sample[,4],probs = c(0.10,0.90))
# #boxplot(output$gibb.sample[,4][c(output$gibb.sample[,4]>qnt[1]) & c(output$gibb.sample[,4] <qnt[2])],main="nu.sq")
# 
# qnt<-quantile(output$gibb.sample[,3],probs = c(0.1,0.9))
# #boxplot(output$gibb.sample[,3][c(output$gibb.sample[,3]>qnt[1]) & c(output$gibb.sample[,3] <qnt[2])], main="sigma.sq")
# 
# qnt<-quantile(output$gibb.sample[,2],probs = c(0.1,0.9))
# #boxplot(output$gibb.sample[,2][c(output$gibb.sample[,2]>qnt[1]) & c(output$gibb.sample[,2] <qnt[2])], main="beta1")
# 
# print(apply(output$gibb.sample,2,median))
# 
# print(apply(output$gibb.sample,2,mean))
# 
