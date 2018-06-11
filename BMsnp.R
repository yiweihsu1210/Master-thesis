#rm(list=ls())
#setwd("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/BMsnp")
 #setwd("~/Dropbox/YiWeiHsu/Rcode/main.code/BMsnp")
# source("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/DV_test.r")
 source("~/Dropbox/YiWeiHsu/Rcode/main.code/DV_test.r")
# source("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/sim.snp.tree.yg.r")
  source("~/Dropbox/YiWeiHsu/Rcode/main.code/sim.snp.tree.yg.r")


LSS.bm<-function(k=k,trait=trait,V=V,D=D){
  n<-length(trait)
  V.inv<-pseudoinverse(V)
  mu.hat<-pseudoinverse(t(D)%*%V.inv%*%D)%*%t(D)%*%V.inv%*%trait
  sigma.sq.hat<- t(trait - D%*%mu.hat)%*%V.inv%*%(trait - D%*%mu.hat)/n
  NegLogML <- (n/2)*log(2*pi)+(1/2)*t(trait- D%*%mu.hat)%*%V.inv%*%(trait- D%*%mu.hat) + (1/2)*log(abs(det(V))) #natural log multivariate noraml PDF???log
  loglik<- -NegLogML
  lss<- 2*loglik-k*log(n)
  return(out=list(lss=lss,mu.hat=mu.hat , sigma.sq.hat = sigma.sq.hat , loglik=loglik))
  }

permute.LSS.bm<- function(k,trait=trait,tree=tree,sims=sims){
  mu.hat.array<-array(0,c(sims))
  sigma.sq.hat.array<-array(0,c(sims))
  loglik.array<-array(0,c(sims))
  lss.array<-array(0,c(sims))
  for(sim.index in 1:sims){
    #    print(sim.index)
    permute.trait<-sample(trait)
    permute.lss.bm <- LSS.bm(k,trait=permute.trait,tree=tree)
    lss.array[sim.index]<-permute.lss.bm$lss
    mu.hat.array[sim.index]<-permute.lss.bm$mu.hat
    sigma.sq.hat.array[sim.index]<-permute.lss.bm$sigma.sq.hat
    loglik.array[sim.index]<-permute.lss.bm$loglik
     }
  return(out= list(lss.array = lss.array , mu.hat.array=mu.hat.array , sigma.sq.hat.array=sigma.sq.hat.array , loglik.array = loglik.array))
  }
