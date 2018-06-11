rm(list=ls())
setwd("/Users/yiweihsu/Documents/BMvsOUSimsHaploid/")
# start.time<-proc.time()
# time.elapse<-proc.time() - start
# print(time.elapse)

filenames<-list.files(path=".")

for (fileIndex in 1:length(filenames)) {
  #fileIndex<-3
  load(filenames[fileIndex])
}

taxa.size.array<-c(10,30,50)
num.tree<-10
num.repdata<-50
p.value.bm.array<-array(0,c(length(taxa.size.array),num.tree,num.repdata))
p.value.ou.array<-array(0,c(length(taxa.size.array),num.tree,num.repdata))

k.bm.array<-array(0,c(length(taxa.size.array),num.tree,num.repdata))
k.ou.array<-array(0,c(length(taxa.size.array),num.tree,num.repdata))

rec.lss.bm.array<-array(0,c(length(taxa.size.array),num.tree,num.repdata))
rec.lss.ou.array<-array(0,c(length(taxa.size.array),num.tree,num.repdata))

sigma.sq.bm.array<-array(0,c(length(taxa.size.array),num.tree,num.repdata))
sigma.sq.ou.array<-array(0,c(length(taxa.size.array),num.tree,num.repdata))

alpha.array<-array(0,c(length(taxa.size.array),num.tree,num.repdata))


for (taxaIndex in 1:length(taxa.size.array)) {
    #taxaIndex<-1
    for (treeIndex in 1:num.tree) {
     #   treeIndex<-2     
        #need to be careful about tree cannot be clustered.
      for (repdataIndex in 1:num.repdata) {
      #   repdataIndex<-1
       # print(paste("taxa" , taxa.size.array[taxaIndex] , "tree" , treeIndex, "data",repdataIndex , sep=""))
      output<-get(paste("taxa" , taxa.size.array[taxaIndex] , "tree" , treeIndex, "data",repdataIndex , sep=""))
      #ls(output)
      if(exists("output") & !is.null(output$bm.DV.data)){
      
      #output$obs.LSS.bm.array
      #output$obs.lss.best.k.ou
      #output$obs.lss.best.k.bm
      
      p.value.bm.array[taxaIndex,treeIndex,repdataIndex]<-output$p.value.bm
      p.value.ou.array[taxaIndex,treeIndex,repdataIndex]<-output$p.value.ou
      
      k.bm.array[taxaIndex,treeIndex,repdataIndex]<-output$best.k.bm
      k.ou.array[taxaIndex,treeIndex,repdataIndex]<-output$best.k.ou
      
      rec.lss.bm.array[taxaIndex,treeIndex,repdataIndex]<-output$obs.lss.best.k.bm$lss
      rec.lss.ou.array[taxaIndex,treeIndex,repdataIndex]<-output$obs.lss.best.k.ou$lss
      
      sigma.sq.bm.array[taxaIndex,treeIndex,repdataIndex]<-output$obs.lss.best.k.bm$sigma.sq.hat
      sigma.sq.ou.array[taxaIndex,treeIndex,repdataIndex]<-output$obs.lss.best.k.ou$sigma.sq.hat
      
      alpha.array[taxaIndex,treeIndex,repdataIndex]<-output$obs.lss.best.k.ou$alpha.hat
      rm(output)
      }#only read non null DV.data analysis
    }
  }
}

taxa10tree2data1
taxa30
p.value.bm.array[1,2,]
p.value.bm.array[2,2,]


#BM and OU are tne same at tree in 1,4,9 not working
par(mfrow=c(2,3))

for(taxaIndex in 1:length(taxa.size.array)){

  print("--------")
  boxplot( c(p.value.bm.array[taxaIndex, !apply(p.value.bm.array[taxaIndex,,], 1, function(y)all(y==0)), ]),main=paste("bm",taxa.size.array[taxaIndex],sep=""))
  boxplot( c(p.value.ou.array[taxaIndex, !apply(p.value.ou.array[taxaIndex,,], 1, function(y)all(y==0)), ]),main=paste("ou",taxa.size.array[taxaIndex],sep=""))
  
  print(summary( c(p.value.bm.array[taxaIndex, !apply(p.value.bm.array[taxaIndex,,], 1, function(y)all(y==0)), ]),main=paste("bm",taxa.size.array[taxaIndex],sep="")))
  print(summary( c(p.value.ou.array[taxaIndex, !apply(p.value.ou.array[taxaIndex,,], 1, function(y)all(y==0)), ]),main=paste("ou",taxa.size.array[taxaIndex],sep="")))
  
  
  }

########################################################
################boxplot for p-value#####################
########################################################

library(ggplot2)
library(reshape2)
ggplot(data = df.m, aes(x=variable, y=value)) + geom_boxplot(aes(fill=Label))


par(mfrow=c(1,3))
for(taxaIndex in 1:length(taxa.size.array)){
A <- c(p.value.bm.array[taxaIndex, !apply(p.value.bm.array[taxaIndex,,], 1, function(y)all(y==0)), ])
B <- c(p.value.ou.array[taxaIndex, !apply(p.value.ou.array[taxaIndex,,], 1, function(y)all(y==0)), ])

mydf <- data.frame(y=c(A,B),x=c(rep("BM",length(A)),rep("OU",length(B))))
with(mydf, boxplot(y~x,main=paste("taxa size",taxa.size.array[taxaIndex],sep="")))
}



df <- rbind(A, B)
df <- rep(c("A", "B"), c(dim(A)[1], dim(B)[1]))
?melt


#####################################################
##############boxplot for sigma.sq###################
#####################################################
library(xtable)
par(mfrow=c(2,3))

for(taxaIndex in 1:length(taxa.size.array)){
  
  print("--------")
  boxplot( c(sigma.sq.bm.array[taxaIndex, !apply(sigma.sq.bm.array[taxaIndex,,], 1, function(y)all(y==0)), ]),main=paste("bm",taxa.size.array[taxaIndex],sep=""))
  boxplot( c(sigma.sq.ou.array[taxaIndex, !apply(sigma.sq.ou.array[taxaIndex,,], 1, function(y)all(y==0)), ]),main=paste("ou",taxa.size.array[taxaIndex],sep=""))
  
  print(summary( c(sigma.sq.bm.array[taxaIndex, !apply(sigma.sq.bm.array[taxaIndex,,], 1, function(y)all(y==0)), ]),main=paste("bm",taxa.size.array[taxaIndex],sep="")))
  print(summary( c(sigma.sq.ou.array[taxaIndex, !apply(sigma.sq.ou.array[taxaIndex,,], 1, function(y)all(y==0)), ]),main=paste("bm",taxa.size.array[taxaIndex],sep="")))
  
  
}

sigma.sq.quantile.array<-array(0,dim = c(2,9))

par(mfrow=c(1,3))
tableIndex<-1
for(taxaIndex in 1:length(taxa.size.array)){
  A <- c(sigma.sq.bm.array[taxaIndex, !apply(sigma.sq.bm.array[taxaIndex,,], 1, function(y)all(y<=0)), ])
  B <- c(sigma.sq.ou.array[taxaIndex, !apply(sigma.sq.ou.array[taxaIndex,,], 1, function(y)all(y<=0)), ])
  
  #print("--------")
  sigma.sq.quantile.array[1,tableIndex:(tableIndex+2)]<-quantile(A, probs = c(0.5,0.025, 0.975))
  sigma.sq.quantile.array[2,tableIndex:(tableIndex+2)]<-quantile(B, probs = c(0.5,0.025, 0.975))
  tableIndex<-tableIndex+3
  
  #mydf <- data.frame(y=c(A,B),x=c(rep("BM",length(A)),rep("OU",length(B))))
  #with(mydf, boxplot(y~x,main=paste("taxa size",taxa.size.array[taxaIndex],sep="")))
}

rownames(sigma.sq.quantile.array)<-c("BM","OU")
colnames(sigma.sq.quantile.array)<-c("10","","","30","","","50","","" )
sigma.sq.quantile.array<-round(sigma.sq.quantile.array,2)
xtable(sigma.sq.quantile.array)


#####################################################
##############boxplot for alpha######################
#####################################################
alpha.quantile.array<-array(0,dim = c(1,9))
tableIndex<-1
for(taxaIndex in 1:length(taxa.size.array)){
  A <-  c(alpha.array[taxaIndex, !apply(alpha.array[taxaIndex,,], 1, function(y)all(y==0)), ])
  
  
  #print("--------")
  alpha.quantile.array[1,tableIndex:(tableIndex+2)]<-quantile(A, probs = c(0.5,0.025, 0.975))
  tableIndex<-tableIndex+3
  
  #mydf <- data.frame(y=c(A,B),x=c(rep("BM",length(A)),rep("OU",length(B))))
  #with(mydf, boxplot(y~x,main=paste("taxa size",taxa.size.array[taxaIndex],sep="")))
}

rownames(alpha.quantile.array)<-"alpha"
colnames(alpha.quantile.array)<-c("10","","","30","","","50","","")
alpha.quantile.array<-round(alpha.quantile.array,2)
xtable(alpha.quantile.array)
