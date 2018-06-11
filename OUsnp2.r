#rm(list=ls())
source("~/Dropbox/YiWeiHsu/Rcode/main.code/DV_test.r")
source("~/Dropbox/YiWeiHsu/Rcode/main.code/sim.snp.tree.yg.r")

library(ape)
library(mvtnorm)
library(invgamma)
library(geiger)
#put OU likelihood
VaF<-function(alpha,V=V){
  A<-exp(-2*alpha*(1-V))
  B<-(1-exp(-2*alpha*V))/(2*alpha)
  return(A*B)
}

# d.mtx<-function(n,snp=NULL){#formula (3) Thompson 2016
#   Z.mtx <- array(0,c(n,2*n))
#   #diplid index
#   for(d.index in 1:n){
#     Z.mtx[d.index,((2*d.index-1):(2*d.index))]<-0.5
#   }
#   return(Z.mtx)
# }

likelihood<-function(model.params,Ye=Ye,Yg=Yg,r=r, k=k,rho=rho,V=V,D=D){
  alpha<-model.params[1]
  beta<-model.params[2:(r+2)]
  sigma.sq<-model.params[(r+3):(r+3)]
  nu.sq<-model.params[(r+4):(r+4)]
  mu<-model.params[(r+4+1):(r+4+k)]
  diploid.size<-length(Ye)
  X<- matrix(runif(n=diploid.size,min=-5,max=5),ncol=1) #shall it be fixed all time ?
  one<-array(1,c(diploid.size,1))
  des.X<-cbind(one,X)
  #Ye.like<-dmvnorm( c(Ye) ,mean= des.X%*%beta, sigma= nu.sq*diag(1,c(length(Ye),length(Ye))),log=TRUE)
  
  #Z<-d.mtx(n=length(Yg))
  Va<-VaF(alpha,V=V)
  #Yg.like <- dmvnorm( c(Yg), mean=Z%*%D%*%mu, sigma=sigma.sq*pseudoinverse(Z%*%Va%*%t(Z)),log=TRUE)
  
  #  Y.loglike<-dmvnorm( c(Yg)+c(Ye), mean = des.X%*%beta + Z%*%D%*%mu ,  sigma = nu.sq*diag(1,c(length(Ye),length(Ye))) + sigma.sq*pseudoinverse(Z%*%Va%*%t(Z)) , log = TRUE)
  Y.loglike<-dmvnorm( c(Yg)+c(Ye), mean = des.X%*%beta + D%*%mu ,  sigma = nu.sq*diag(1,c(length(Ye),length(Ye))) + sigma.sq*pseudoinverse(Va), log = TRUE)
  return( Y.loglike )
}


prior<-function(model.params,r=r,k=k,rho=rho){
  alpha<-model.params[1]
  beta<-model.params[2:(r+2)]
  sigma.sq<-model.params[(r+3):(r+3)]
  nu.sq<-model.params[(r+4):(r+4)]
  mu<-model.params[(r+4+1):(r+4+k)]

  alpha.prior<-dexp(alpha,rate=1,log=TRUE)

  u.sq<-100
  u.sq.I<-u.sq*diag(1,c(r+1,r+1))
  beta0<-array(0,c(r+1,1))
  beta.prior<-dmvnorm(beta, mean=beta0, sigma=u.sq.I,log=TRUE) #PUT MVNNORM HERE

  #rho<-0.5 #will do a function input later
  a<-6
  b<-200*rho^2
  if(b==0){
    sigma.sq.prior <- dinvgamma(sigma.sq,shape=a,scale=1e-5,log=TRUE)
    }else{
    sigma.sq.prior <- dinvgamma(sigma.sq,shape=a,scale=b,log=TRUE)}

  c<-6
  d<-200*(1-rho)^2
  if(d==0){
    nu.sq.prior<-dinvgamma(nu.sq,shape=c,scale=1e-5,log=TRUE)
    }else{
    nu.sq.prior<-dinvgamma(nu.sq,shape=c,scale=d,log=TRUE)
    }
  mu0<-array(0,c(k,1))
  w.sq<-100
  w.sq.I<-w.sq*diag(1,c(k,k))
  mu.prior<-dmvnorm(mu,mean=mu0,sigma=w.sq.I,log=TRUE)

  return(alpha.prior+beta.prior+sigma.sq.prior+nu.sq.prior+mu.prior)
}


posterior<-function(model.params,Yg=Yg,Ye=Ye,r=r,k=k,rho=rho,V=V,D=D){
  return(likelihood(model.params,Yg=Yg,Ye=Ye,r=r,k=k,rho=rho,V=V,D=D) + prior(model.params,r=r,k=k,rho=rho))
}


alpha.proposal<-function(alpha){
  # we can use MLE estimate and sd from geiger packages
  #?fitContinuous
  #trait<-matrix(trait,ncol=1)
  #tree<- rcoal(length(trait))#It would be better if we use the tree we construct. #compute.brlen(stree(length(trait),type="left"))
  #rownames(trait)<-tree$tip.label
  #ou.geiger<-fitContinuous(phy=tree,dat=trait, model="OU")
  #tree<-stree(length(trait),type="star")#can do this, need to fix DV function
  #ou.result<-optimize(OUnegloglike, lower=0,upper=10  ,k=k, trait=trait, tree=tree)
  #alpha
  #use star tree
  #return(rnorm(1,mean= ou.geiger$opt$alpha, sd=ou.geiger$opt$alpha/3 ) )
  return(rnorm(1, mean=alpha, sd = alpha/4))
}



Invgamma.setup<-function(t.mean=t.mean,t.var=t.var){
  shape<- 2+t.mean^2/t.var
  scale<-t.mean*(t.mean^2/t.var+1)
  return(list(shape=shape,scale=scale))
}



MHGibbs<-function(Yg=Yg,Ye=Ye,k=k,V=V,D=D,X=X, rho=rho,r=r, true.beta=true.beta,true.sigma.sq = true.sigma.sq,nu.sq=true.nu.sq , sim.chains=sim.chains,a=a,b=b,c=c,d=d,e=e){
  # sim.chains<-10
  # rho<-0.5
  # diploid.size<-5
  # k=3
  Y=Yg+Ye
  #diploid.size<-length(Y)
  r<-1
  u.sq<- 0.001#20#100
  u.sq.I<-u.sq*diag(1,c(r+1,r+1))
  #print(true.beta)
  beta0<-array(c(true.beta[1],true.beta[2]),c(r+1,1))
  beta<-matrix(rmvnorm(n=1,mean=beta0,sigma=u.sq.I),ncol=1)#(14)
  sigma.sq<-1#(16)
  # a<-6
  # b<-0.2 #200*rho^2
  # c<-6
  # d<-10  #200*(1-rho)^2
  #nu.sq<-rinvgamma(1,shape=c,scale=d)
#  I<-diag(1,c(diploid.size))
  I<-diag(1,c(taxa.size))
  w.sq<-100
  w.sq.I<- w.sq*diag(1,c(k,k))
  mu0<- array(0,c(k,1))
  mu<-matrix(rmvnorm(n=1,mean=mu0,sigma=w.sq.I),ncol=1)

  alpha<-rexp(1,rate=1/e)

  
  
  #nu.sq=2
  
  
  
  model.params<-c(alpha,beta, sigma.sq, nu.sq, mu )
  names(model.params)<-c( "alpha",paste("b", 0:(length(beta)-1),sep="" ), "sig.sq","nu.sq",paste("mu",1:length(mu),sep=""))
  mhgibb.sample<-array(0,c(sim.chains,length(model.params)))
  colnames(mhgibb.sample)<-c( "alpha",paste("b", 0:(length(beta)-1),sep="" ), "sig.sq","nu.sq",paste("mu",1:length(mu),sep=""))
  mhgibb.sample[1,]<-model.params
  Va<- VaF(alpha,V=V)

  #Z<-d.mtx(n=length(Y))
# ZVatZ.inv<-pseudoinverse(Z%*%Va%*%t(Z))
  Vat.inv<-pseudoinverse(Va)
  one<-array(1,c(length(Y),1))
  des.X<-cbind(one,X)
  tdes.X.des.X<-t(des.X)%*%des.X
  #ZD<-Z%*%D



  for(i in 2:sim.chains){
    #print(i)
    if(i%%100==0)print(i)
    ###alpha<-chain[i-1,1] #need to check

    Y <- Yg+Ye

    #Gibb starts
    sigma.sq.shape <- a + taxa.size
    sigma.sq.scale <- b +(0.5*t(Yg-D%*%mu)%*%Vat.inv%*%(Yg-D%*%mu))
    sigma.sq <- rinvgamma(1,shape=sigma.sq.shape, scale = 1/sigma.sq.scale)
    #sigma.sq<-rnorm(1,mean=true.sigma.sq, sd= sqrt( sigma.sq.scale^2/ ( (sigma.sq.shape-1 )^2*(sigma.sq.shape-2))))
   
    
    nu.sq.shape<-c + taxa.size
    nu.sq.scale<-d + (0.5*t(Y- des.X%*%beta-Yg)%*%(Y- des.X%*%beta-Yg))
    nu.sq<-rinvgamma(1,shape=nu.sq.shape,scale=1/nu.sq.scale)
    #nu.sq<-rnorm(1,mean=true.nu.sq, sd= sqrt(nu.sq.scale^2/((nu.sq.shape-1)^2*(nu.sq.shape-2))))
    
    # V.yg<-pseudoinverse(ZVatZ.inv/sigma.sq  + diag(1,c(diploid.size,diploid.size))/nu.sq)
    # V.yg.stuff<-V.yg%*%((Y-des.X%*%beta)/nu.sq  + ZVatZ.inv%*%ZD%*%mu/sigma.sq)
    # Yg<-matrix(rmvnorm(n=1,mean=V.yg.stuff,sigma=V.yg),ncol=1)
    
    
    #V.yg<-pseudoinverse(Vat.inv/sigma.sq  + diag(1,c(diploid.size,diploid.size))/nu.sq)
    V.yg<-pseudoinverse(Vat.inv/sigma.sq  + diag(1,c(taxa.size,taxa.size))/nu.sq)     
    V.yg.stuff<-V.yg%*%((Y-des.X%*%beta)/nu.sq  + Vat.inv%*%D%*%mu/sigma.sq)
    Yg<-matrix(rmvnorm(n=1,mean=V.yg.stuff,sigma=V.yg),ncol=1)
    
    
    # V.mu <- pseudoinverse( t(ZD)%*%ZVatZ.inv%*%ZD/sigma.sq + diag(1,c(k,k))/w.sq )
    # V.mu.stuff<-V.mu%*%(  t(ZD)%*%ZVatZ.inv%*%Yg /sigma.sq + mu0/w.sq  )
    # mu<-matrix(rmvnorm(1,mean=V.mu.stuff,sigma=V.mu),ncol=1)
    V.mu <- pseudoinverse( t(D)%*%Vat.inv%*%D/sigma.sq + diag(1,c(k,k))/w.sq )
    V.mu.stuff<-V.mu%*%(  t(D)%*%Vat.inv%*%Yg /sigma.sq + mu0/w.sq  )
    mu<-matrix(rmvnorm(1,mean=V.mu.stuff,sigma=V.mu),ncol=1)
    
    u.sq<- 0.001
    V.beta<-pseudoinverse(tdes.X.des.X/nu.sq + diag(1,dim(tdes.X.des.X))/u.sq )
    V.beta.stuff<-V.beta%*%( t(des.X)%*%(Y-Yg) /nu.sq + beta0/u.sq )
    # print(V.beta)
    # print(eigen(V.beta)$values)
    beta<-matrix(rmvnorm(n=1,mean=V.beta.stuff, sigma=V.beta), ncol=1)

    
    
    #i=1
    mhgibb.sample[i,2:(length(model.params))]<-c(beta, sigma.sq,nu.sq,mu)
    #print(paste("mhgibb.sample is  ",mhgibb.sample[i,2:(length(model.params))],sep=""))
    #Gibb ends
    #mh starts
    alpha.pps<-alpha.proposal(mhgibb.sample[i-1,1]) #slower
    #print(paste("alpha.pps is  ",alpha.pps,sep=""))
    
    proposal<-c(alpha.pps, mhgibb.sample[i,2:(length(model.params))])
    
    #print(proposal)
    #print(paste("propasal is  ",proposal,sep=""))
    #    likelihood(proposal,Yg=Yg,Ye=Ye,r=r,k=k,rho=rho,tree=tree)
    #    model.params<-proposal


    probab=exp(posterior(proposal,Yg=Yg,Ye=Ye,r=r,k=k, rho=rho,V=V,D=D) - posterior(mhgibb.sample[i-1,],Yg=Yg,Ye=Ye,r=r,k=k, rho=rho,V=V,D=D))
    #print(posterior(proposal,Yg=Yg,Ye=Ye,r=r,k=k, rho=rho,V=V,D=D))
    #print(posterior(mhgibb.sample[i-1,],Yg=Yg,Ye=Ye,r=r,k=k, rho=rho,V=V,D=D))
    
    #mhgibb.sample[i,]<-proposal
    
    #print(probab)
    if(runif(1)<probab){
      mhgibb.sample[i,]=c(alpha.pps, mhgibb.sample[i,2:length(model.params)])
    }else{
      mhgibb.sample[i,]=c(mhgibb.sample[i-1,1], mhgibb.sample[i,2:length(model.params)] )
    }
    #mh ends

    X<- matrix(runif(n=taxa.size,min=-10,max=10),ncol=1)
    Ze<-matrix(rmvnorm(n=1,mean=beta[1]+beta[2]*X,sigma=sqrt(abs(sigma.sq))*diag(1,c(taxa.size, taxa.size))),ncol=1)# may use differnt a,b,s for  a + bX, s*diag...
    Ye<-(1-rho)*Ze
  }#end for loop
  return(list(mhgibb.sample=mhgibb.sample))
}
