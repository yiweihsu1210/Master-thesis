rm(list=ls())
setwd("/Users/yiweihsu/Documents/BMsnp2/")
#load("simtest1.RData")
# start.time<-proc.time()
# time.elapse<-proc.time() - start
# print(time.elapse)

#filenames<-list.files(path=".")

# for (fileIndex in 1:length(filenames)) {
#   load(filenames[fileIndex])
# }

# taxa.size.array<-c(10,30,50)
# taxa.size.array<-c(50)
# for (taxaIndex in 1:length(taxa.size.array)) {
#   print(paste("BMdiploidSimsTaxa",taxa.size.array[taxaIndex],".RData" , sep=""))
#   load.Rdata(paste("BMdiploidSimsTaxa",taxa.size.array[taxaIndex],".RData" , sep=""))
#   print(taxa.size.array)
#   taxa.size.array<-c(10,30,50)
#   }
#setwd("/Users/yiweihsu/Documents/BMsnp2/")
#save.image("simtest1.RData")


ls(taxa20tree3data10rho0.2k1)
head(taxa20tree3data10rho0.2k1$gibb.sample)
ls()

taxa.size.array<-c(20,30,50)
num.tree<-10
num.repdata<-50
rho.array<-seq(0,1,0.2)
k.array<- 3:5
sim.chains<-5000#50000


output.arrayk1<-array(0,c(length(taxa.size.array),num.tree,num.repdata,length(rho.array),sim.chains,7))
output.arrayk2<-array(0,c(length(taxa.size.array),num.tree,num.repdata,length(rho.array),sim.chains,8))
output.arrayk3<-array(0,c(length(taxa.size.array),num.tree,num.repdata,length(rho.array),sim.chains,9))

existed.count<-0
total.count<-0
for (taxaIndex in 1:length(taxa.size.array)) {
  #taxaIndex<-1
#  load(paste("BMdiploidSimsTaxa",taxa.size.array[taxaIndex],".RData" , sep=""))
#  taxa.size.array<-c(20,30,50)
  taxa.size<-taxa.size.array[taxaIndex]
  for (treeIndex in 1:num.tree) {
     #treeIndex<-1
      for (repdataIndex in 1:num.repdata) {
      repdataIndex<-1
        for(rhoIndex in 1:length(rho.array)){
        #rhoIndex<-1
          for (kIndex in 1:length(k.array)) {
            #kIndex<-1
            total.count<-total.count+1
          print(paste("taxa" , taxa.size , "tree" , treeIndex, "data",repdataIndex,"rho",  rho.array[rhoIndex], "k", kIndex,sep=""))
          
          
          if(kIndex==1){
            object2get<-paste("taxa" , taxa.size , "tree" , treeIndex, "data",repdataIndex,"rho",  rho.array[rhoIndex], "k", kIndex,sep="")
            if(exists(object2get)){
              existed.count<-existed.count+1
              try(output <- get(object2get))
              if(!is.null(output)){
                output.arrayk1[taxaIndex,treeIndex,repdataIndex,rhoIndex,,]<-output$gibb.sample
                print(output$gibb.sample[1,])
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
                output.arrayk2[taxaIndex,treeIndex,repdataIndex,rhoIndex,,]<-output$gibb.sample
                print(output$gibb.sample[1,])
              }
            }
          }
          if(kIndex==3){
            #gibb.sample<-array(0,c(sim.chains,9))
            #output<-list(gibb.sample=gibb.sample)
            object2get<-paste("taxa" , taxa.size , "tree" , treeIndex, "data",repdataIndex,"rho", rho.array[rhoIndex], "k", kIndex,sep="")
            if(exists(object2get)){
              existed.count<-existed.count+1
              try(output <- get(object2get))
              if(!is.null(output)){
                output.arrayk3[taxaIndex,treeIndex,repdataIndex,rhoIndex,,]<-output$gibb.sample
                print(output$gibb.sample[1,])
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

setwd("/Users/yiweihsu/Documents/BMsnp2/")
#save (taxa.size.array , num.tree, num.repdata , rho.array  , k.array  , sim.chains , output.arrayk1 , output.arrayk2 , output.arrayk3, file = "summaryoutput.RData")
###################################################
##############SummaryPart II#######################
###################################################
rm(list=ls())
setwd("/Users/yiweihsu/Documents/BMsnp2/")
load("summaryoutput.RData")

#dim(output.arrayk1)
#taxa,tree,data,rho,sim,param
#3   10   50    6 5000    7
#param
#beta0, beta1, sigma.sq, nu.sq, mu1, mu2,mu3

taxaIndex<-1  # 10, 30, 50
treeIndex<-3  # 1:10
dataIndex<-10  # 1:50
rhoIndex<-2  # 0, 0.2,0.4,0.6,0.8,1 
paramIndex<-2 # beta0, beta1, sigma.sq, nu.sq, mu1, mu2,mu3

boxplot(output.arrayk1[taxaIndex,treeIndex,dataIndex,rhoIndex,,paramIndex])
summary(output.arrayk1[taxaIndex,treeIndex,dataIndex,rhoIndex,,paramIndex])
quantile(output.arrayk1[taxaIndex,treeIndex,dataIndex,rhoIndex,,paramIndex], prob=seq(0,1,by=0.1))




sim.chains.red<-sim.chains/100
k1rho<-array(0,c( length(taxa.size.array), num.tree, num.repdata*sim.chains.red, length(rho.array)))

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
x <- k1rho_test
typeof(x)
head(x)
x$"0"

y <- sapply(1:6,function(i){
  remove_outliers(x[,i])
})

