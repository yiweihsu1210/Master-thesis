#setwd("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/OUsnp")
#setwd("~/Dropbox/YiWeiHsu/Rcode/main.code/OUsnp")

#rm(list=ls())
#source("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/DV_test.r")
source("~/Dropbox/YiWeiHsu/Rcode/main.code/DV_test.r")
#source("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/sim.snp.tree.yg.r")
source("~/Dropbox/YiWeiHsu/Rcode/main.code/sim.snp.tree.yg.r")
VaF<-function(alpha,V=V){
  A<-exp(-2*alpha*(1-V))
  B<- (1-exp(-2*alpha*V))/(2*alpha)
  return(A*B)
  }

#for OU we need optimization here
OUnegloglike<-function(alpha, k=k, trait=trait, V=V,D=D){
    Va<-VaF(alpha,V=V)
    Va.inv<-pseudoinverse(Va)
    mu.hat<-pseudoinverse(t(D)%*%Va.inv%*%D)%*%t(D)%*%Va.inv%*%trait
    sigma.sq.hat<- t(trait - D%*%mu.hat)%*%Va.inv%*%(trait - D%*%mu.hat)/length(trait)
#    print(mu.hat)
#    print(sigma.sq.hat)
    NegLogML <- (length(trait)/2)*log(2*pi)+(1/2)*t(trait- D%*%mu.hat)%*%Va.inv%*%(trait- D%*%mu.hat) + (1/2)*log(abs(det(Va)))
    return(NegLogML)
    }

MLEsOU <- function(alpha, k=k, trait=trait, V=V, D=D){
  Va<-VaF(alpha,V=V)
  Va.inv<-pseudoinverse(Va)
  mu.hat<-pseudoinverse(t(D)%*%Va.inv%*%D)%*%t(D)%*%Va.inv%*%trait
  sigma.sq.hat<- t(trait - D%*%mu.hat)%*%Va.inv%*%(trait - D%*%mu.hat)/length(trait)
  return(list(mu.hat=mu.hat,sigma.sq.hat=sigma.sq.hat,alpha.hat=alpha))
  }

LSS.ou<-function(k,trait=trait,V=V, D=D){
  n<-length(trait)
  ou.result<-optimize(OUnegloglike, lower=0.01,upper=10 , k=k, trait=trait, V=V, D=D)
  #print(ou.result)
  #ou.mle<-MLEsOU(c(ou.result$objective), k=k, trait=trait, tree=tree)
  loglik<- -ou.result$objective
  lss<- 2*loglik - k*log(n)
  alpha.hat = ou.result$minimum
  mlesou<-MLEsOU(alpha = alpha.hat , k=k , trait = trait, V=V , D=D)
  return(out=list(lss=lss , loglik=loglik , alpha.hat = alpha.hat, mu.hat=mlesou$mu.hat,sigma.sq.hat=mlesou$sigma.sq.hat))
  }

permute.LSS.ou<- function(k,trait=trait,tree=tree,sims=sims){
  lss.array<-array(0,c(sims))
  for(sim.index in 1:sims){
    print(sim.index)
    permute.trait<-sample(trait)
    permute.lss.ou <- LSS.ou(k,trait=permute.trait,tree=tree)
    lss.array[sim.index]<-permute.lss.ou
    }
  return(lss.array)
  }
