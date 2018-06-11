#move back to maincode folder after implementation
rm(list=ls())
source("~/Dropbox/YiWeiHsu/Rcode/main.code/OUsnp2_haploid/OUsnp2.r")
setwd("~/Dropbox/YiWeiHsu/Rcode/simulation/OUsnpGEhaploid/Taxa50/")

rho.array<- c(0.1,0.3,0.5,0.7,1)#,0.25,0.3,0.35,0.5,0.7,0.9,1)
true.b1.array<- c(4,3.9,3.8,3.7,3.6)#,3,2.8,2.7,2.1,1.3,0.7,0.1)
true.b0.array<-c(10,10,10,10,10)#rep(10,length(rho.array))
true.sigma.sq.array<-c(20,22,24,26,28)# seq(20,100,10)
true.nu.sq.array<- c(250,255,260,265,270)#seq(250,410,20)
#true.sigma.sq.array<-true.nu.sq.array

true.alpha.array<- c(0.01,0.15,0.2,0.5,1)#c(0.1,0.5,1,1.5,2,2.5,3,3.5,5)
length(true.alpha.array)

#diploid.size.array <- c(10)
taxa.size.array <- c(50)#2*diploid.size.array
k.array<- c(5)#3:4 #3:5
r<-1
seqlen<-1000
# ouparams<-c(0.5,0.5,5,5,90,80,100)
# #ouparams<-c(0.001,0.001,200,200,90,80,100)
# names(ouparams)<-c("alp1","alp2","sig1","sig2","th0","th1","th2")
num.tree<-1#10
num.repdata<-1#50
sim.chains<-100#50000
# true.sigma.sq <- 5
# true.b0=20
# true.b1=40
# true.nu.sq=2
# true.beta<-c(true.b0,true.b1)

for (taxaIndex in 1:length(taxa.size.array)) {
  #taxaIndex <-1
  taxa.size<-taxa.size.array[taxaIndex]
  #diploid.size<-taxa.size/2
  #generate a tree 
  #taxa.size<-taxa.size.array[1]
  for (treeIndex in 1:num.tree) {
    #treeIndex<-1
    tree<-sim.snp.tree(taxa.size=taxa.size,seqlen=seqlen)
    tree<-chronos(tree)
    #plot(tree)
    tree$tip.label<-paste("t",tree$tip.label,sep="")
    
    #tree<-pbtree(n=taxa.size)
    
    vcv.tree<-vcv(tree)
    diag(vcv.tree)<-max(vcv.tree)
    vcv.tree<-vcv.tree/max(vcv.tree)
    diag(vcv.tree)<-diag(vcv.tree)
    try(tree<-vcv2phylo(vcv.tree))
    plot(tree)
    
    for(repdataIndex in 1:num.repdata){
      #repdataIndex<-1
      for(rhoIndex in 1:length(rho.array)){
        #rhoIndex<-1
        print("using different rho")
        rho<-rho.array[rhoIndex]
        true.b0<-true.b0.array[rhoIndex]
        true.b1<-true.b1.array[rhoIndex]
        true.sigma.sq <- true.sigma.sq.array[rhoIndex]
        true.nu.sq<- true.nu.sq.array[rhoIndex]
        true.beta<-c(true.b0,true.b1)
        true.alpha<- true.alpha.array[rhoIndex]
        ouparams<-c(true.alpha,true.alpha,true.sigma.sq,true.sigma.sq,90,80,100)
        names(ouparams)<-c("alp1","alp2","sig1","sig2","th0","th1","th2")
        #generate trait : (1) genetic component (2) environmental component
        Zg.X <- sim.snp.yg(tree=tree,ouparams=ouparams)#trait only
        #Z<-d.mtx(n=diploid.size.array[1])
        Zg<-matrix(Zg.X$X,ncol=1)
        #rownames(Zg) <- Zg.X$Genus_species
        Yg<-rho*Zg
        #print(dim(Yg))
        #print(paste("iid estimate for sigma2 is  ",var(Yg),sep=""))
        
        #problem where nu.sq variance so bigm due to the multivaraite normal centered on different position
        X<- matrix(runif(n=taxa.size,min=-5,max=5),ncol=1)
        #print(X)
        #print(true.nu.sq*diag(1,c(taxa.size, taxa.size)))
        Ze<-matrix(rmvnorm(n=1,mean=true.b0+true.b1*X,sigma=true.nu.sq*diag(1,c(taxa.size, taxa.size))),ncol=1)# may use differnt a,b,s for  a + bX, s*diag...
        #print(paste("true.nu.sq is ",  true.nu.sq,"  ,iid estimate for nu.sq is  ", var(Ze),sep=""))
        
        #Ze<-matrix(rmvnorm(n=1,mean= array(200,c(taxa.size,1)),sigma=true.nu.sq*diag(1,c(taxa.size, taxa.size))),ncol=1)# may use differnt a,b,s for  a + bX, s*diag...
        
        #print(dim(Ze))
        #print(paste("true.nu.sq is ",  true.nu.sq,"  ,iid estimate for nu.sq is  ", var(Ze),sep=""))
        
        
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
            
            ab<-Invgamma.setup(t.mean=true.sigma.sq,t.var=1)
            a<-ab$shape-taxa.size
            a
            b<-ab$scale
            cd<-Invgamma.setup(t.mean=true.nu.sq,t.var=1)
            c<-cd$shape-taxa.size
            d<-cd$scale
            
            e<-true.alpha
            
            try(output<-MHGibbs(Yg=Yg,Ye=Ye, k=k.array[kIndex],V=V, D=D,X=X, rho=rho,r=r, true.beta=true.beta,true.sigma.sq = true.sigma.sq,nu.sq=true.nu.sq ,sim.chains=sim.chains,a=a,b=b,c=c,d=d,e=e))
            
            assign(paste("taxa" , taxa.size , "tree" , treeIndex, "data",repdataIndex ,"rho",rho,"k",k.array[kIndex], sep=""),output)
            #boxplot(output$gibb.sample[,2], main= "beta1")
            #boxplot(output$gibb.sample[,3], main="sigma.sq")
            #boxplot(output$gibb.sample[,4], main="nu.sq")
            
          }#end of DVisNULL
        }#end of kIndex
      }#end of rhoIndex
      setwd("/Users/yiweihsu/Documents/OUsnp2v2/")
      #save.image(file="OUGEhaploidSimsTaxa50_1.RData")
      setwd("~/Dropbox/YiWeiHsu/Rcode/simulation/OUsnpGEhaploid/Taxa50/")
    }#end of repdataIndex
  }#end of treeIndex
}#end of taxaIndex

alpha.rho.array<-array(0,c(sim.chains,length(rho.array)))
colnames(alpha.rho.array)<-rho.array
b0.rho.array<-array(0,c(sim.chains,length(rho.array)))
colnames(b0.rho.array)<-rho.array
b1.rho.array<-array(0,c(sim.chains,length(rho.array)))
colnames(b1.rho.array)<-rho.array
sigma.sq.rho.array<-array(0,c(sim.chains,length(rho.array)))
colnames(sigma.sq.rho.array)<-rho.array
nu.sq.rho.array<-array(0,c(sim.chains,length(rho.array)))
colnames(nu.sq.rho.array)<-rho.array
sigma.nu.rho.array<-array(0,c(sim.chains,length(rho.array)))
colnames(sigma.nu.rho.array)<-rho.array

#output.array<-array(0,c(length(taxa.size.array),num.tree,num.repdata,length(rho.array),sim.chains,9))
for (taxaIndex in 1:length(taxa.size.array)) {
  taxa.size<-taxa.size.array[taxaIndex]
  #taxaIndex <-1
  for (treeIndex in 1:num.tree) {
    #treeIndex<-1
    for(repdataIndex in 1:num.repdata){
      #repdataIndex<-1
      for(rhoIndex in 1:length(rho.array)){
        rho<-rho.array[rhoIndex]
        #rhoIndex<-1
        for(kIndex in 1:length(k.array)){
          #kIndex<-1
          #assign(paste("taxa" , taxa.size , "tree" , treeIndex, "data",repdataIndex,"rho", rho,"k", k.array[kIndex], sep=""),output)
          output<-get(paste("taxa" , taxa.size , "tree" , treeIndex, "data",repdataIndex,"rho", rho, "k", k.array[kIndex],sep=""))  
          alpha.rho.array[,rhoIndex]<-output$mhgibb.sample[,1]
          b0.rho.array[,rhoIndex]<-output$mhgibb.sample[,2]
          b1.rho.array[,rhoIndex]<-output$mhgibb.sample[,3]
          sigma.sq.rho.array[,rhoIndex]<-output$mhgibb.sample[,4]
          nu.sq.rho.array[,rhoIndex]<-output$mhgibb.sample[,5]
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

setwd("/Users/yiweihsu/Dropbox/YiWeiHsu/翊瑋")
jpeg(file="ouboxplot.jpeg")
par(mfrow=c(2,3))
plot(tree)
boxplot(apply(alpha.rho.array,2,remove_outliers),xlab=expression(rho),ylab=expression(alpha))
points(1:5,true.alpha.array,col=2,pch=13)
boxplot(apply(b1.rho.array,2,remove_outliers),xlab=expression(rho),ylab=expression(beta[1]))
points(1:5,true.b1.array,col=2,pch=13)
boxplot(apply(sigma.sq.rho.array,2,remove_outliers),,xlab=expression(rho),ylab=expression(sigma^2))
points(1:5,true.sigma.sq.array,col=2,pch=13)
boxplot(remove_outliers(nu.sq.rho.array),,xlab=expression(rho),ylab=expression(nu^2))
points(1:5,true.nu.sq.array,col=2,pch=13)
boxplot(sigma.nu.rho.array,,xlab=expression(rho),ylab=expression(sigma^2/(sigma^2+nu^2)))
points(1:5,true.sigma.sq.array/(true.sigma.sq.array+true.nu.sq.array),col=2,pch=13)
dev.off()
#return posterior from gibbs sampling
#par(mfrow=c(1,3))
#qnt<-quantile(output$gibb.sample[,4],probs = c(0.10,0.90))
#boxplot(output$gibb.sample[,4][c(output$gibb.sample[,4]>qnt[1]) & c(output$gibb.sample[,4] <qnt[2])],main="nu.sq")

#qnt<-quantile(output$gibb.sample[,3],probs = c(0.1,0.9))
#boxplot(output$gibb.sample[,3][c(output$gibb.sample[,3]>qnt[1]) & c(output$gibb.sample[,3] <qnt[2])], main="sigma.sq")

#qnt<-quantile(output$gibb.sample[,2],probs = c(0.1,0.9))
#boxplot(output$gibb.sample[,2][c(output$gibb.sample[,2]>qnt[1]) & c(output$gibb.sample[,2] <qnt[2])], main="beta1")
#print(apply(output$mhgibb.sample,2,median))
