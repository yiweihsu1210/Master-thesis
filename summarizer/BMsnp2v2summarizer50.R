rm(list=ls())
setwd("/Users/yiweihsu/Documents/BMsnp2v2/")
load("BMGEhaploidSimsTaxa50_1.RData")

# rho.array<- c(0.0001,0.1,0.25,0.3,0.35,0.5,0.7,0.9,1)
# taxa.size.array <- c(50)# 2*diploid.size.array
# k.array<- c(5)
# r<-1
# seqlen<-1000
# num.tree<-1#10
# num.repdata<-1#50
# sim.chains<-50000


b0.rho.array<-array(0,c(sim.chains,9))
colnames(b0.rho.array)<-rho.array
b1.rho.array<-array(0,c(sim.chains,9))
colnames(b1.rho.array)<-rho.array
sigma.sq.rho.array<-array(0,c(sim.chains,9))
colnames(sigma.sq.rho.array)<-rho.array
nu.sq.rho.array<-array(0,c(sim.chains,9))
colnames(nu.sq.rho.array)<-rho.array
sigma.nu.rho.array<-array(0,c(sim.chains,9))
colnames(sigma.nu.rho.array)<-rho.array

#output.array<-array(0,c(length(taxa.size.array),num.tree,num.repdata,length(rho.array),sim.chains,9))
for (taxaIndex in 1:length(taxa.size.array)) {
  #taxaIndex <-1
  for (treeIndex in 1:num.tree) {
    #treeIndex<-1
    for(repdataIndex in 1:num.repdata){
      #repdataIndex<-1
      for(rhoIndex in 1:length(rho.array)){
        #rhoIndex<-1
        for(kIndex in 1:length(k.array)){
          #kIndex<-1
          output<-get(paste("taxa" , taxa.size , "tree" , treeIndex, "data",repdataIndex,"rho", rhoIndex, "k", k.array[kIndex],sep=""))  
          b0.rho.array[,rhoIndex]<-output$gibb.sample[,1]
          b1.rho.array[,rhoIndex]<-output$gibb.sample[,2]
          sigma.sq.rho.array[,rhoIndex]<-output$gibb.sample[,3]
          nu.sq.rho.array[,rhoIndex]<-output$gibb.sample[,4]
          sigma.nu.rho.array[,rhoIndex]<-sigma.sq.rho.array[,rhoIndex]/(sigma.sq.rho.array[,rhoIndex]+nu.sq.rho.array[,rhoIndex])
        }
      }
    }
  }
}

remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75),na.rm = na.rm)
  H <- 1.5 * IQR(x)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
  # if(y<-NA){
  #   try(y)
  # }
}

par(mfrow=c(2,2))
boxplot(apply(b1.rho.array,2,remove_outliers),xlab=expression(rho),ylab=expression(beta[1]))
boxplot(apply(sigma.sq.rho.array,2,remove_outliers),,xlab=expression(rho),ylab=expression(sigma^2))
boxplot(remove_outliers(nu.sq.rho.array),,xlab=expression(rho),ylab=expression(nu^2))
boxplot(sigma.nu.rho.array,,xlab=expression(rho),ylab=expression(sigma^2/(sigma^2)))
