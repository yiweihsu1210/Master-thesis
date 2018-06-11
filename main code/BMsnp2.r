#move back to maincode folder after implementation
# source("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/DV_test.r")
# source("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/YiWeiHsu/Rcode/main.code/sim.snp.tree.yg.r")
source("~/Dropbox/YiWeiHsu/Rcode/main.code/DV_test.r")
source("~/Dropbox/YiWeiHsu/Rcode/main.code/sim.snp.tree.yg.r")

library(corpcor)
library(ape)
library(mvtnorm)
library(invgamma)

# d.mtx<-function(n,snp=NULL){ #formula (3) Thompson 2016
#   Z.mtx<-array(0,c(n,2*n))
#   #diploid index
#   for(d.index in 1:n){
#     Z.mtx[d.index,((2*d.index-1):(2*d.index))]<-0.5
#     Z.mtx
#   }
#   return(Z.mtx)
# }

Invgamma.setup<-function(t.mean=t.mean,t.var=t.var){
  shape<- 2+t.mean^2/t.var
  scale<-t.mean*(t.mean^2/t.var+1)
  return(list(shape=shape,scale=scale))
  }

Gibbs<-function(Yg=Yg,Ye=Ye, k=k,V=V, D=D,X=X, rho=rho, r=r, true.beta = true.beta, sigma.sq = sigma.sq, nu.sq=nu.sq,sim.chains=sim.chains,a=a,b=b,c=c,d=d){
  
  
  Y=Yg+Ye

  #diploid.size<-length(Y)
  r<-1
  u.sq<- 0.001
  u.sq.I<-u.sq*diag(1,c(r+1,r+1))
  beta0<-array(c(true.beta[1],true.beta[2]),c(r+1,1))
  beta<-matrix(rmvnorm(n=1,mean=beta0,sigma=u.sq.I),ncol=1)#(14)
  
  #a<-6
  #b<-0.2#1/(a-1) #200*rho^2 = 200*0.0001=0.02
  #sigma.sq<-rinvgamma(1,shape = a,scale = b)
  
  #c<-6
  #d<-10 #1/(200*(1-rho)^2) #2*(c-1)=d
  #nu.sq<-rinvgamma(1,shape=c,scale=d)
  #  nu.sq<-rinvgamma(1,shape=c,rate=d)
  #print("first samle from prior")
  #print(c(sigma.sq,nu.sq))
  
  #I<-diag(1,c(diploid.size))
  I<-diag(1,c(taxa.size))
  w.sq<-1
  w.sq.I<- w.sq*diag(1,c(k,k))
  
  mu0<- array(0,c(k,1))
  mu<-matrix(rmvnorm(n=1,mean=mu0,sigma=w.sq.I),ncol=1)
  
  model.params<-c(beta, sigma.sq, nu.sq, mu)
  #  print(model.params)
  names(model.params)<-c( paste("b", 0:(length(beta)-1),sep="" ), "sig.sq","nu.sq",paste("mu",1:length(mu),sep=""))
  gibb.sample<-array(0,c(sim.chains,length(model.params)))
  colnames(gibb.sample)<-c( paste("b", 0:(length(beta)-1),sep="" ), "sig.sq","nu.sq",paste("mu",1:length(mu),sep=""))
  gibb.sample[1,]<-model.params
  
  Vt.inv<-pseudoinverse(V)
  
  one<-array(1,c(length(Y),1))
  des.X<-cbind(one,X)
  tdes.X.des.X<- t(des.X)%*%des.X
  #ZD<-Z%*%D
  
  
  for(chainIndex in 2:sim.chains){
    #chainIndex<-2
#    print("Y")
#    print(Y)
    Y <- Yg+Ye
    if(chainIndex %%100 ==0){print(chainIndex)}
    
    sigma.sq.shape <-     a + taxa.size
  
    #b<-200*rho^2
    sigma.sq.scale <-  b +(0.5*t(Yg-D%*%mu)%*%Vt.inv%*%(Yg-D%*%mu))
#    print("var term")
#    print(t(Yg-D%*%mu)%*%Vt.inv%*%(Yg-D%*%mu))   
    #    print(mu)
#    print(cbind(Yg,D%*%mu))
    
#    print( paste("second term", 0.5*t(Yg-mean(Yg))%*%(Yg-mean(Yg)), sep=" "))
    #print( paste("second term", 0.5*t(Yg-D%*%mu)%*%Vt.inv%*%(Yg-D%*%mu), sep=" "))
#    print(paste("here is sigma.scale",sigma.sq.scale,sep=" "))
#    sigma.sq<-rinvgamma(1, shape=sigma.sq.shape,scale=1/72)#(19)  #sig.sq | y,yg,beta,mu,nu.sq
    
    sigma.sq<-rinvgamma(1, shape=sigma.sq.shape,scale=1/sigma.sq.scale)#(19)  #sig.sq | y,yg,beta,mu,nu.sq
    #sigma.sq<-rinvgamma(1,shape=18,scale=1/72)
    
    
    #sigma.sq<-true.sigma.sq
    
    #c=7
    nu.sq.shape <- c + taxa.size
    #d=200*(1-rho)^2
    
    nu.sq.scale <- d + (0.5*t(Y- des.X%*%beta-Yg)%*%(Y- des.X%*%beta-Yg))
#    print(paste("here is nu.sq.scale",nu.sq.scale,sep=" "))
          
#    nu.sq<-rinvgamma(1,shape=nu.sq.shape,scale=1/72) #(20)
    nu.sq<-rinvgamma(1,shape=nu.sq.shape,scale=1/nu.sq.scale) #(20)
    #nu.sq<-rinvgamma(1,shape=18,scale=1/72)
    
    #nu.sq<-true.nu.sq
    
    
    
    
    V.yg<-pseudoinverse(Vt.inv/sigma.sq + diag(1,c( taxa.size,taxa.size))/nu.sq)#(22)
    V.yg.stuff<-V.yg%*%( (Y-des.X%*%beta)/nu.sq  +  Vt.inv%*%D%*%mu/sigma.sq)
    Yg<-matrix(rmvnorm(n=1,mean=V.yg.stuff,sigma=V.yg),ncol=1) #(21)
#    print("Yg")
#    print(Yg)
    
    V.mu <- pseudoinverse(t(D)%*%Vt.inv%*%D/sigma.sq + diag(1,c(k,k))/w.sq )#(24)
    V.mu.stuff <-  V.mu%*% (  t(D)%*%Vt.inv%*%Yg  /sigma.sq   +    mu0/w.sq)
    mu<- matrix(rmvnorm(n=1,mean= V.mu.stuff, sigma=V.mu),ncol=1)#(23)
    
    u.sq<- 0.001
    V.beta<-pseudoinverse(  tdes.X.des.X /nu.sq + diag(1,dim(tdes.X.des.X) )/u.sq)#(26)
    V.beta.stuff<-V.beta%*%( t(des.X)%*%(Y-Yg)/nu.sq +  beta0 /u.sq )
    #print("beta center")
    #print(V.beta.stuff)
    #print(paste("betamean   ", round(V.beta.stuff,2),sep=""))
    beta <- matrix(rmvnorm(n=1, mean = V.beta.stuff, sigma=V.beta),ncol=1)    #(25)
    #print(cbind(V.beta.stuff,beta0))
    #print(V.beta)
    
    #beta<-matrix(rmvnorm(n=1,mean=beta0,sigma=u.sq.I),ncol=1)#(14)
    
    #beta<-beta0
    
    X<- matrix(runif(n=taxa.size,min=-10,max=10),ncol=1)
    Ze<-matrix(rmvnorm(n=1,mean=beta[1]+beta[2]*X,sigma=sqrt(abs(sigma.sq))*diag(1,c(taxa.size, taxa.size))),ncol=1)# may use differnt a,b,s for  a + bX, s*diag...
    Ye<-(1-rho)*Ze
    
 #   print("Ye")
#    print(Ye)
    
    gibb.sample[chainIndex,]<-c(beta, sigma.sq, nu.sq, c(mu))
  }#end for loop
  #print(gibb.sample[,2])
  #plot(gibb.sample[,2])
  return(list(gibb.sample=gibb.sample))
}#end of function
#return posterior from gibbs sampling
