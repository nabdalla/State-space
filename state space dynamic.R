#Bayesian dynamic model
#general model: measurement equation: Y(t)=C(t)+vt vt iid Pv
#transition equation: C(t+1)=(1-Q'+K.lV/V)C(t)+G'/V+wt, wt iid Pw
#Q'=Q+erf QR+Ql, G'=(1-el)G
#Simulation
library(nimble)
library(rjags)
library(plotrix)
library(compositions)
library(mvtnorm)
library(expm)
library(DPpackage)
n=100
Qprime=13.8
Gprime=351.54
V=100
sigma<-0.5
aomega<-2
bomega<-1
kl=0.1
set.seed(2017)
wt<-rgamma(n,aomega,bomega)

#el<-3
#Gprime<-(1-el)*G
yt<-NULL
c<-NULL
#values of c very huge or very small and hence yt so taking log is problem
for(i in 1:(n-1)){
  c[1]<-1
  c[i+1]<-(1-(Qprime+kl*V)/V)*c[i]+Gprime/V+wt[i+1]
}
plot(c)
Qprime/V+kl
set.seed(000)
vt<-rnorm(n,0,sigma)
#logyt<-log(c)+(vt)
#yt<-exp(logyt)
#hist(log(yt))
logyt<-log(c)+(vt)
plot((logyt))
hist(exp(logyt))
# p1<-p<-yhat<-k<-NULL
# for(i in 1:(n-1)){
# p1[1]<-1
# p[i]<-((1-(Qprime+kl*V)/V)^2)*p1[i]+(aomega/bomega^2)
# k[i]<-(p[i]*(1/c[i]))/((1/c[i])^2*p[i]+sigma)
# p1[i+1]<-p[i]-k[i]*p[i]/c[i]
# yhat[i]<-c[i]+k[i]*(exp(logyt[i])-1)}
# plot(yhat)
# points(exp(logyt), col="red")
# points((c), col="blue")
#boostrap draft (particle filter)
#create function for pf
pf<-function(nparticles,logyt, Qprime, Gprime, kl,aomega,bomega,sigma){
  #initialize vectors
  cestsir<-cpartiold<-cpartinew<-wparti<-wtotal<-wpartin<-cpartiw<-cpartires<-wpartires<-
    cestpf<-stdsir<-phat<-phat1<-NULL
  nparti<-nparticles
  nparti1<-1/nparti
  #distribution initial time
  cestsir[1]<-logyt[1]
  for(i in 1:nparti){
    cpartiold[i]<-exp(logyt[1])+rgamma(1,aomega,bomega)
  }
  #advance particles
  for(j in 1:(n-1)){
    for(i in 1:nparti){
      #prediction
      cpartinew[i]<-(1-(Qprime+kl)/V)*cpartiold[i]+rgamma(1,aomega,bomega)
      #weights using likelihood
      wparti[i]<-exp(-((log(cpartinew[i])-logyt[j+1])/sigma)^2)
    } 
    wtotal<-sum(wparti)
    wpartin<-wparti/wtotal
    phat<-(wtotal/nparti)
    phat1<-c(phat1,phat)
    
    #update
    cpartiw<-cpartinew*wpartin
    #mean at time j+1
    #cestpf[j+1]<-cpartinew*wparti
    #resample
    cresa<-rep(0,nparti)
    uresa<-rep(0,nparti)
    #cumulative sum of weights
    cresa[1]<-wpartin[1]
    for(i in 2:nparti){
      cresa[i]<-cresa[i-1]+wpartin[i]
    }
    iresa<-1
    uresa[1]<-runif(1,0,1)*nparti1
    for(k in 1:nparti){
      uresa[k]<-uresa[1]+nparti1*(k-1)
      while (uresa[k]>cresa[iresa]){
        iresa=iresa+1}
      cpartires[k]<-cpartinew[iresa]
      wpartires[k]<-nparti1
    }
    cestsir[j+1]<-mean(cpartires)
    stdsir[j+1]<-sqrt(var(cpartires))
    cpartiold<-cpartires
    #end of resampling
    #cpartiold<-cpartinew
    mat<-(cestsir)
  }
  updates<-as.numeric(mat)
  return(c(updates,phat1))}
truepf<-pf(200,logyt, Qprime, Gprime, kl,aomega,bomega,sigma)
cpf<-truepf[1:100]
likpf<-truepf[101:200]
plot(cpf)
points(exp(logyt), col="red")
points(c,col="blue")
legend("bottomright",col=c("black","red","blue"),text.font=4,legend=
         c("Particle filter", "Yt","Ct"),bg="white", lty=1,cex=0.5)
#mcmc
n.iters<-1000
tune<-0.01
thin<-10
theta<-c(Qprime, Gprime, kl, aomega,bomega,sigma)
p<-length(theta)
ll<- -1e99
theta.mat<-matrix(0,nrow=n.iters, ncol=p)
x.mat<-matrix(0,nrow=n.iters,ncol=100)
pfmcmc1<-function(n.iters, tune, thin,theta){
  for(i in 1:n.iters){
    message(paste(i,""),appendLF=FALSE)
    for(j in 1:thin){
      thetaprop<-theta*exp(rnorm(p,0,tune))
      priorprop<-log(c(dunif(thetaprop[1],11,17),dunif(thetaprop[2],281,482),dunif(thetaprop[3],0,0.8),
                       dunif(thetaprop[4],0.5,3),dunif(thetaprop[5],0.5,3),dgamma(thetaprop[6],1,2)))
      prior<-log(c(dunif(theta[1],11,17),dunif(theta[2],281,482),dunif(theta[3],0,0.8),
                   dunif(theta[4],0.5,3),dunif(theta[5],0.5,3),dgamma(theta[6],1,2)))
      pfprop<-pf(200,logyt,thetaprop[1],thetaprop[2],thetaprop[3],thetaprop[4],thetaprop[5],thetaprop[6])
      llprop<-log(prod((pfprop)[101:200],na.rm=TRUE))
      pftheta<-pf(200,logyt,theta[1],theta[2],theta[3],theta[4],theta[5],theta[6])
      x<-pftheta[1:100]
      ll<-log(prod((pftheta)[101:200],na.rm=TRUE))
      xprop<-pfprop[1:100]
      if(log(runif(1))< llprop+priorprop[1]+priorprop[2]+priorprop[3]+
         priorprop[4]+priorprop[5]+priorprop[6]-(ll+prior[1]+prior[2]+prior[3]+
                                                 prior[4]+prior[5]+prior[6])){
        theta<-thetaprop
        ll=llprop
        x<-xprop[1:100]
      } else{
        theta<-theta
        ll=ll
        x<-x
      }
    }
    theta.mat[i,]<-theta
    x.mat[i,]<-x
  }
  return(list(theta.mat,x.mat))}
ptm <- proc.time()
pfmcmc.out<-pfmcmc1(n.iters,tune, thin,theta)
proc.time() - ptm

thetajulia<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmc.csv")
apply(thetajulia, 2, function(x){quantile(x)})
xjulia<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc.csv")
pred.julia<-apply(xjulia, 2, mean)

#library(compiler)
#pc<-cmpfun(pfmcmc1)
#ptm <- proc.time()
#f<-pc(n.iters,tune, thin,theta)
#proc.time() - ptm

apply(pfmcmc.out[[1]], 2, function(x){quantile(x)})
quartz()
mcmcSummary(pfmcmc.out[[1]])
pred<-apply(pfmcmc.out[[2]], 2, mean)
quartz()
plot(pred,c)
abline(0,1)
quartz()
plot(cpf)
#points(pred, col="red")
points(c,col="blue")
points(exp(logyt),col="green")
points(pred.julia, col="red")
legend("bottomright",col=c("black","red","blue","green"),text.font=4,legend=
         c("Particle filter","PMMH","Ct","yt"),bg="white", lty=1,cex=0.8)
quartz()
acf(pfmcmc.out[[1]][,4], lag.max=1000)
#jags
cat("model{
    for (i in 1:N){ 
    omega[i+1]~dgamma(aomega,bomega)
    c[i+1]<-(1-(Qprime+kl*V)/V)*c[i]+Gprime/V+omega[i+1]
    logyt[i+1]~dnorm(log(c[i+1]), sigma)
    # p[i]<-((1-(Qprime+kl*V)/V)^2)*p1[i]+(aomega/bomega^2)
    # k[i]<-p[i]/(p[i]+prec)
    # p1[i+1]<-p[i]-k[i]*p[i]
    #logyhat[i+1]~dnorm(k[i+1]*log(c[i])+(k[i+1]*c[i])-(k[i+1]*c[i]),p[i+1])
    }
    #p1[1]<-1
    c[1]<-1
    logyt[1]~dnorm(log(c[1]), sigma)
    sigma<-1/prec
    aomega~dunif(1,3)
    bomega~dunif(0.5,1.5)
    prec~dgamma(2,1)
    Qprime~dunif(11,17)
    Gprime~dunif(281,482)
    kl~dunif(0,0.8)
    }",file="modeljags.jag")
jags<-jags.model("modeljags.jag",data=list("logyt"=logyt[1:100], N=99, "V"=V))
update(jags, 1000)
mcmcjags<-jags.samples(jags,
                       c('Qprime','Gprime',"kl",'sigma','c','aomega','bomega'),
                       1000)
mcmcsamplesjags<-coda.samples(jags,
                              c('Qprime','Gprime','kl','sigma','c','aomega','bomega'),
                              1000)
quartz()
plot(mcmcsamplesjags)
summary(mcmcsamplesjags)
m1.mcmc<-(as.mcmc(mcmcsamplesjags))
m1.mat<-as.matrix(mcmcsamplesjags)
m1.dat<-as.data.frame(m1.mat)
aomega.post1zone<-m1.dat$aomega
bomega.post1zone<-m1.dat$bomega
qprime.post1zone<-m1.dat$Qprime
gprime.post1zone<-m1.dat$Gprime
kl.post1zone<-m1.dat$kl
sigma.post1zone<-1/m1.dat$sigma
quantile(sigma.post1zone)
c.hat <- apply(mcmcjags$c, 1, mean)
c.hatlower<- apply(mcmcjags$c, 1, function(x) quantile(x,0.025))
c.hatupper<-apply(mcmcjags$c, 1, function(x) quantile(x,0.975))
quartz()
plotCI(c, c.hat,ui=c.hatupper, li=c.hatlower)
abline(c(0,1))
pf.post1zone<-vector("list")
for(i in 1:1000){
  pf.post1zone[[i]]<-pf(200,logyt,qprime.post1zone[i],gprime.post1zone[i],kl.post1zone[i],aomega.post1zone[i],bomega.post1zone[i],sigma.post1zone[i])
}
logyhat<-matrix(0,nrow=100,1000)
for(i in 1:1000){
  for(j in 1:100){
    logyhat[j,i]<-(pf.post1zone[[i]][j])
  }
}
logyhat.1zone<-apply(logyhat,1,mean)
logyt.hatlower<- apply(logyhat, 1, function(x) quantile(x,0.025))
logyt.hatupper<-apply(logyhat, 1, function(x) quantile(x,0.975))
quartz()
plotCI(truepf,logyhat.1zone, ui=logyt.hatupper, li=logyt.hatlower)
abline(c(0,1))
quartz()
#plot(truepf, col="red",ylim=c(0,30))
#points((logyhat.1zone))
plot(c,col="red")
points(c.hatupper,col="blue")
points(c.hatlower,col="blue")
points(c.hat, col="purple")
points(exp(logyt),col="green")
legend("bottomright",col=c("red","blue","blue","purple","green"),text.font=4,legend=
         c("Ct","Chatupper","Chatlower","chat","yt"),bg="white", lty=1,cex=0.8)
plot(truepf,logyhat.1zone)
abline(0,1)
#DP

pf.dp1<-function(nparticles,logyt, Qprime, Gprime, kl,aomega,bomega){
  #initialize vectors
  cestsir<-cpartiold<-cpartinew<-wparti<-wtotal<-wpartin<-cpartiw<-cpartires<-wpartires<-
    cestpf<-stdsir<-NULL
  nparti<-nparticles
  nparti1<-1/nparti
  #distribution initial time
  cestsir[1]<-logyt[1]
  for(i in 1:nparti){
    cpartiold[i]<-exp(logyt[1])+rgamma(1,aomega,bomega)
  }
  #advance particles
  for(j in 1:(n-1)){
    for(i in 1:nparti){
      #prediction
      cpartinew[i]<-(1-(Qprime+kl)/V)*cpartiold[i]+rgamma(1,aomega,bomega)
      #weights using likelihood
      for(l in 1:100){
        g01<-rnorm(100,0,1)
        g02<-rgamma(100,2,1)
        b1<-rbeta(100,1,1)
        p1<-numeric(100)
        p1[1]<-b1[1]
        p1[2:100]<-sapply(2:100, function(i) b1[i] * prod(1 - b1[1:(i-1)]))
        mu<-sample(g01,prob=p1,size=nparticles,replace=TRUE)
        b2<-rbeta(100,1,1)
        p2<-numeric(100)
        p2[1]<-b2[1]
        p2[2:100]<-sapply(2:100, function(i) b2[i] * prod(1 - b2[1:(i-1)]))
        sigma<-sample(g02,prob=p2,size=nparticles,replace=TRUE)
      }
      wparti[i]<-exp(-0.5*(logyt[j+1]-log(cpartinew[i])-mu[i])^2/sigma[i])
    } 
    wtotal<-sum(wparti)
    wpartin<-wparti/wtotal
    #update
    cpartiw<-cpartinew*wpartin
    #mean at time j+1
    #cestpf[j+1]<-cpartinew*wparti
    #resample
    cresa<-rep(0,nparti)
    uresa<-rep(0,nparti)
    #cumulative sum of weights
    cresa[1]<-wpartin[1]
    for(i in 2:nparti){
      cresa[i]<-cresa[i-1]+wpartin[i]
    }
    iresa<-1
    uresa[1]<-runif(1,0,1)*nparti1
    for(k in 1:nparti){
      uresa[k]<-uresa[1]+nparti1*(k-1)
      while (uresa[k]>cresa[iresa]){
        iresa=iresa+1}
      cpartires[k]<-cpartinew[iresa]
      wpartires[k]<-nparti1
    }
    cestsir[j+1]<-mean(cpartires)
    stdsir[j+1]<-sqrt(var(cpartires))
    cpartiold<-cpartires
    #end of resampling
    #cpartiold<-cpartinew
    mat<-(cestsir)
  }
  updates<-as.numeric(mat)
  return(updates)}
dp.true<-pf.dp1(200,logyt, Qprime, Gprime, kl,aomega,bomega)
plot(c)
points(updates,col="red")
points(exp(logyt),col="blue")
#jags
cat("model{
    for (i in 1:N){  
    omega[i+1]~dgamma(aomega,bomega)
    c[i+1]<-(1-(Qprime+kl*V)/V)*c[i]+Gprime/V+omega[i+1]
    logyt[i+1]~dnorm(log(c[i+1])+muv[z1[i+1]],sigmav[z2[i+1]])
    z1[i+1]~dcat(p1[])
    z2[i+1]~dcat(p2[])
    #pt[i]<-((1-(Qprime+kl*V)/V)^2)*p1t[i]+(aomega/bomega^2)
    # #Kalman filter gain
    # k[i]<-pt[i]/(pt[i]+sigma)
    # #covariance
    # p1t[i+1]<-pt[i]-k[i]*pt[i]
    # #update
    # logyhat[i]<-k[i]*logyt[i]
    }
    #p1t[1]<-1
    p1[1]<-r1[1]
    p2[1]<-r2[1]
    for(j in 2:(m-1)){p1[j]<-r1[j]*(1-r1[j-1])*p1[j-1]/r1[j-1]
    p2[j]<-r2[j]*(1-r2[j-1])*p2[j-1]/r2[j-1]}
    for(l in 1:(m-1)){r1[l]~dbeta(1,alpha)
    r2[l]~dbeta(1,alpha)}
    ps1<-sum(p1[1:(m-1)])
    p1[m]<-1-ps1
    ps2<-sum(p2[1:(m-1)])
    p2[m]<-1-ps2
    for(l in 1:m){
    muv[l]~dnorm(basemu, basetau)
    sigmav[l]~dgamma(2,1)
    }
    basemu<-0
    basetau<-1
    c[1]<-1
    logyt[1]~dnorm(log(c[1]), prec)
    alpha<-1
    aomega~dunif(1,3)
    bomega~dunif(0.5,1.5)
    sigma~dgamma(2,1)
    prec<-1/sigma
    Qprime~dunif(11,17)
    Gprime~dunif(281,482)
    kl~dunif(0,0.8)
    }",file="modeljagsdp.jag")
jags.dp<-jags.model("modeljagsdp.jag",data=list("logyt"=logyt[1:41], N=40, "V"=V,"m"=100))
update(jags.dp, 1000)
mcmcjags.dp<-jags.samples(jags.dp,
                          c('Qprime','Gprime','kl','sigma','c','aomega','bomega'),
                          1000)
mcmcsamplesjags.dp<-coda.samples(jags.dp,
                                 c('Qprime','Gprime','kl','sigma','c','aomega','bomega'),
                                 1000)
quartz()
plot(mcmcsamplesjags.dp)
summary(mcmcsamplesjags.dp)
m1.mcmc.dp<-(as.mcmc(mcmcsamplesjags.dp))
m1.mat.dp<-as.matrix(mcmcsamplesjags.dp)
m1.dat.dp<-as.data.frame(m1.mat.dp)
aomega.post1zone.dp<-m1.dat.dp$aomega
bomega.post1zone.dp<-m1.dat.dp$bomega
qprime.post1zone.dp<-m1.dat.dp$Qprime
gprime.post1zone.dp<-m1.dat.dp$Gprime
kl.post1zone.dp<-m1.dat.dp$kl

c.hat.dp <- apply(mcmcjags.dp$c, 1, mean)
c.hatlower.dp<- apply(mcmcjags.dp$c, 1, function(x) quantile(x,0.025))
c.hatupper.dp<-apply(mcmcjags.dp$c, 1, function(x) quantile(x,0.975))
quartz()
plotCI(c[1:41], c.hat.dp,ui=c.hatupper.dp, li=c.hatlower.dp)
abline(c(0,1))
pf.post1zone.dp<-vector("list")
for(i in 1:1000){
  pf.post1zone.dp[[i]]<-pf.dp1(100,logyt,qprime.post1zone.dp[i],gprime.post1zone.dp[i],kl.post1zone.dp[i],aomega.post1zone.dp[i],bomega.post1zone.dp[i])
}
logyhat.dp<-matrix(0,nrow=41,3)
for(i in 1:3){
  for(j in 1:41){
    logyhat.dp[j,i]<-(pf.post1zone.dp[[i]][j])
  }
}
logyhat.1zone.dp<-apply(logyhat.dp,1,mean)
logyt.hatlower.dp<- apply(logyhat.dp, 1, function(x) quantile(x,0.025))
logyt.hatupper.dp<-apply(logyhat.dp, 1, function(x) quantile(x,0.975))
quartz()
plotCI(dp.true[1:41],logyhat.1zone.dp, ui=logyt.hatupper.dp, li=logyt.hatlower.dp)
abline(c(0,1))
quartz()
plot(dp.true[1:41], col="red",ylim=c(0,30))
points((logyhat.1zone.dp))
points(c[1:41],col="blue")
points(c.hat.dp, col="red")
points(exp(logyt[1:41]),col="green")

plot(dp.true[1:41],logyhat.1zone.dp)
abline(0,1)

#DP all
#jags
cat("model{
    for (i in 1:N){ 
    omega[i+1]~dgamma(aomega[z1[i+1]],bomega[z2[i+1]])
    c[i+1]<-(1-(Qprime+kl*V)/V)*c[i]+Gprime/V+omega[i+1]
    logyt[i+1]~dnorm(log(c[i+1])+muv[z3[i+1]],sigmav[z4[i+1]])
    z1[i+1]~dcat(p1[])
    z2[i+1]~dcat(p2[])
    z3[i+1]~dcat(p3[])
    z4[i+1]~dcat(p4[])    
    pt[i]<-((1-(Qprime+kl*V)/V)^2)*p1t[i]+(aomega[z1[i+1]]/bomega[z2[i+1]]^2)
    k[i]<-pt[i]/(pt[i]+(1/sigmav[z4[i+1]]))
    p1t[i+1]<-pt[i]-k[i]*pt[i]
    logyhat[i]<-k[i]*logyt[i]
    }
    p1t[1]<-1
    p1[1]<-r1[1]
    p2[1]<-r2[1]
    p3[1]<-r3[1]
    p4[1]<-r4[1]
    for(j in 2:(m-1)){p1[j]<-r1[j]*(1-r1[j-1])*p1[j-1]/r1[j-1]
    p2[j]<-r2[j]*(1-r2[j-1])*p2[j-1]/r2[j-1]
    p3[j]<-r3[j]*(1-r3[j-1])*p3[j-1]/r3[j-1]
    p4[j]<-r4[j]*(1-r4[j-1])*p4[j-1]/r4[j-1]}
    for(l in 1:(m-1)){r1[l]~dbeta(1,alpha)
    r2[l]~dbeta(1,alpha)
    r3[l]~dbeta(1,alpha)
    r4[l]~dbeta(1,alpha)
    }
    ps1<-sum(p1[1:(m-1)])
    p1[m]<-1-ps1
    ps2<-sum(p2[1:(m-1)])
    p2[m]<-1-ps2
    ps3<-sum(p3[1:(m-1)])
    p3[m]<-1-ps3
    ps4<-sum(p4[1:(m-1)])
    p4[m]<-1-ps4
    for(l in 1:m){
    muv[l]~dnorm(basemu, basetau)
    sigmav[l]~dgamma(2,1)
    aomega[l]~dunif(1,5)
    bomega[l]~dunif(0.5,2)
    }
    basemu<-0
    basetau<-1
    c[1]<-1
    logyt[1]~dnorm(log(c[1]), sigmav[1])
    alpha<-1
    Qprime~dunif(11,17)
    Gprime~dunif(281,482)
    kl~dunif(0,0.8)
    }",file="modeljagsdp2.jag")

jags.dp2<-jags.model("modeljagsdp2.jag",data=list("logyt"=logyt[1:41], N=40, "V"=V,"m"=100))
update(jags.dp2, 1000)
mcmcjags.dp2<-jags.samples(jags.dp2,
                           c('Qprime','Gprime','k','kl','sigmav','c','aomega','bomega','logyhat'),
                           1000)
mcmcsamplesjags.dp2<-coda.samples(jags.dp2,
                                  c('Qprime','Gprime','k','kl','sigmav','c','aomega','bomega','logyhat'),
                                  1000)
quartz()
plot(mcmcsamplesjags.dp2)
summary(mcmcsamplesjags.dp2)

c.hat.dp2 <- apply(mcmcjags.dp2$c, 1, mean)
c.qdp2<-apply(mcmcjags.dp2$c, 1, quantile)
c.hatlowerdp2<- apply(mcmcjags.dp2$c, 1, function(x) quantile(x,0.025))
c.hatupperdp2<-apply(mcmcjags.dp2$c, 1, function(x) quantile(x,0.975))
quartz()
plotCI(c[1:41], c.hat.dp2,ui=c.hatupperdp2, li=c.hatlowerdp2)
abline(0,1)
quartz()
plot((c[1:41]))
points(c.hat.dp2,col="red")
points(c.hat, col="blue")
points(c.hat.dp,col="green")
logy.hatdp2 <- apply(mcmcjags.dp2$logyhat, 1, mean)
logyt.hatlowerdp2<- apply(mcmcjags.dp2$logyhat, 1, function(x) quantile(x,0.025))
logyt.hatupperdp2<-apply(mcmcjags.dp2$logyhat, 1, function(x) quantile(x,0.975))
quartz()
plot((logy.hatdp2))

quartz()

plotCI(logyhat[1:40], logy.hatdp2,ui=logyt.hatupper, li=logyt.hatlower)
abline(c(0,1))

plot(exp(logy.hatdp2),col="blue")
points(exp(logyhat),col="green")

#two compartment model simulation
n=100
Qprime=13.8
Gprime=351.54
VF=1
VN=1
Sigma<-diag(0.1,2)
kl=0.1
beta=5
V=diag(1,2)
r=3
set.seed(2017)
wt<-rlnorm.rplus(n,c(0,0),diag(2))
h=0.01
A<-matrix(c(-(beta)/VN, (beta)/VN, beta/VF, -(beta+Qprime)/VF+kl),nrow=2, ncol=2,
          byrow=TRUE)
eigen(A)
g<-matrix(c(Gprime/VN,0),nrow=1, ncol=2)

yt<-NULL
c<-matrix(c(rep(0,2*n)), nrow=n, ncol=2)
#values of c very huge or very small and hence yt so taking log is problem
for(i in 1:(n-1)){
  c[1,]<-c(0,0.5)
  c[i+1,]<-(c[i,])%*%(h*A+diag(1,2))+(g*h)+h*wt[i+1,]
}
c
plot(c[,1]/VN,ylab="mg/m3", xlab="time")
points(c[,2])

#exact 
c2<-matrix(c(rep(0,2*n)), nrow=n, ncol=2)
#imp

for(i in 1:(n-1)){
  c2[1,]<-c(0,0.5)
  c2[i+1,]<- expm(h*i*A)%*%c2[1,]+solve(A)%*%(expm(h*i*A)-diag(2))%*%t(g)+wt[i+1,]
}
plot(c2[,1])
points(c2[,2])
plot(c[,1],c2[,1])
abline(0,1)
plot(c[,2],c2[,2])
abline(0,1)
(Gprime/beta+Gprime/Qprime)
(Gprime/Qprime)
c<-c2
#or
# c3<-matrix(c(rep(0,2*n)), nrow=n, ncol=2)
# lambda1<-(eigen(A)$value[1])
# lambda2<-(eigen(A)$value[2])
# e.vec11<-eigen(A)$vectors[1,1]
# e.vec12<-eigen(A)$vectors[2,1]
# e.vec21<-eigen(A)$vectors[1,2]
# e.vec22<-eigen(A)$vectors[2,2]
# c22<- 1/(e.vec12-(e.vec11*(e.vec12-e.vec22)/(e.vec11-e.vec12)))
# c11<- -((e.vec12-e.vec22)/((e.vec11-e.vec12)))*c22
# c4<-matrix(c(rep(0,2*n)), nrow=n, ncol=2)
# #imp
# for(i in 1:n){
#   c4[i,]<-c11*(e.vec1)*exp(lambda1*i)+c22*(e.vec2)*exp(lambda2*i)
# }
# for(i in 1:n){
#   c3[i,1]<-(Gprime/Qprime)+(Gprime/beta)+Gprime*(((beta*Qprime)+(lambda2*VN*(beta+Qprime)))/
#                                               (beta*Qprime*VN*(lambda1-lambda2)))*
#     exp(lambda1*i)+Gprime*(((beta*Qprime)+(lambda1*VN*(beta+Qprime)))/
#                                (beta*Qprime*VN*(lambda1-lambda2)))*
#     exp(lambda2*i)
#   c3[i,2]<-Gprime/Qprime+Gprime*((lambda1*VN)+beta)/beta*((beta*Qprime+lambda2*VN*(beta+Qprime))/
#                                                (beta*Qprime*VN*(lambda1-lambda2)))*
#     exp(lambda1*i)+Gprime*(((lambda2*VN)+beta)/beta)*(((beta*Qprime)+lambda1*VN*(beta+Qprime))/
#                              (beta*Qprime*VN*(lambda1-lambda2)))*
#     exp(lambda2*i)
#   
# }
# 
# plot(c2[,1],c3[,1])
# abline(0,1)
# plot(c2[,2],c3[,2])
# abline(0,1)
set.seed(123)
vt<-rmvnorm(n,c(0,0),Sigma)
logyt<-log(c)+(vt)
exp(logyt)
#KF
# p1<-p<-logyhat.2zone<-k<-vector("list")
# for(i in 1:(n-1)){
#   p1[[1]]<-diag(2)
#   p[[i]]<-A%*%p1[[i]]%*%t(A)+diag(2)
#   k[[i]]<-p[[i]]%*%solve(p[[i]]+diag(2))
#   p1[[i+1]]<-p[[i]]-k[[i]]%*%p[[i]]
#   logyhat.2zone[[i]]<-k[[i]]%*%logyt[i,]}
# logyhat1.2zone<-logyhat2.2zone<-NULL
# for(i in 1:99){
#   logyhat1.2zone[i]<-logyhat.2zone[[i]][1,]
#   logyhat2.2zone[i]<-logyhat.2zone[[i]][2,]
# }
# plot((logyhat2.2zone))
#create function for pf
pf2zone<-function(nparticles,n,logyt, A,g,aomega,bomega,Sigma){
  #initialize vectors
  cestsir<-matrix(0,n,2)
  cpartiold<-matrix(0,nparticles,2)
  wparti<-cpartinew<-cpartires<-wpartires<-matrix(0,nparticles,2)
  wtotal1<-wtotal2<-wpartin1<-wpartin2<-cpartiw<-
    cestpf<-stdsir1<-stdsir2<-NULL
  nparti<-nparticles
  nparti1<-1/nparti
  #distribution initial time
  cestsir[1,]<-exp(logyt[1,])
  for(i in 1:nparti){
    cpartiold[i,]<-exp(logyt[1,])+rgamma(2,aomega,bomega)
  }
  #advance particles
  for(j in 1:(n-1)){
    for(i in 1:nparti){
      #prediction
      cpartinew[i,]<-cpartiold[i,]%*%(h*A+diag(1,2))+(g*h)+h*rgamma(2,aomega,bomega)
      #weights using likelihood
      wparti[i,]<-exp(-(t(log(cpartinew[i,])-logyt[j+1,])%*%solve(Sigma)%*%(log(cpartinew[i,])-logyt[j+1,]))/2)
    } 
    wtotal1<-sum(wparti[,1])
    wtotal2<-sum(wparti[,2])
    wpartin1<-wparti[,1]/wtotal1
    wpartin2<-wparti[,2]/wtotal2
    #update
    cpartiw1<-cpartinew[,1]*wpartin1
    cpartiw2<-cpartinew[,2]*wpartin2
    #mean at time j+1
    #cestpf[j+1]<-cpartinew*wparti
    #resample
    cresa1<-rep(0,nparti)
    uresa1<-rep(0,nparti)
    cresa2<-rep(0,nparti)
    uresa2<-rep(0,nparti)
    #cumulative sum of weights
    cresa1[1]<-wpartin1[1]
    cresa2[1]<-wpartin2[1]
    for(i in 2:nparti){
      cresa1[i]<-cresa1[i-1]+wpartin1[i]
      cresa2[i]<-cresa2[i-1]+wpartin2[i]
    }
    iresa1<-1
    iresa2<-1
    uresa1[1]<-runif(1,0,1)*nparti1
    uresa2[1]<-runif(1,0,1)*nparti1
    for(k in 1:nparti){
      uresa1[k]<-uresa1[1]+nparti1*(k-1)
      uresa2[k]<-uresa2[1]+nparti1*(k-1)
      while (uresa1[k]>cresa1[iresa1]){
        iresa1=iresa1+1}
      cpartires[k,1]<-cpartinew[iresa1,1]
      wpartires[k,1]<-nparti1
      while (uresa2[k]>cresa2[iresa2]){
        iresa2=iresa2+1}
      cpartires[k,2]<-cpartinew[iresa2,2]
      wpartires[k,2]<-nparti1
    }
    cestsir[j+1,1]<-mean(cpartires[,1])
    cestsir[j+1,2]<-mean(cpartires[,2])
    stdsir1[j+1]<-sqrt(var(cpartires[,1]))
    stdsir2[j+1]<-sqrt(var(cpartires[,2]))
    cpartiold<-cpartires
    #end of resampling
    #cpartiold<-cpartinew
    mat1<-(cestsir[,1])
    mat2<-(cestsir[,2])
  }
  updates1<-as.numeric(mat1)
  updates2<-as.numeric(mat2)
  return(c(updates1,updates2))
}
ptm <- proc.time()

truepf.2zone<-pf2zone(200,100,logyt, A,g,aomega,bomega,Sigma)
proc.time()-ptm
quartz()
plot(truepf.2zone[1:100])
points(exp(logyt[,1]), col="red")
points((c[,1]), col="blue")
legend("bottomright",col=c("black","red","blue"),text.font=4,legend=
         c("Particle filter true N","measurement yt","CtN"),bg="white", lty=1,cex=0.8)
#jags
quartz ()
plot(truepf.2zone[101:200],col="purple")
points(exp(logyt[,2]), col="red")
points((c[,2]), col="blue")
legend("bottomright",col=c("purple","red","blue"),text.font=4,legend=
         c("Particle filter true F","measurement yt","CtF"),bg="white", lty=1,cex=0.8)
#jags
I<-diag(2)

cat("model{
    for (i in 1:N){
    logomega[i+1,1:2]~dmnorm(c(0,0),vari)
    c[i+1,1]<-c[i,1]*((a*h)+1)+c[i,2]*cc*h+g1*h+h*exp(logomega[i+1,1])
    c[i+1,2]<-c[i,1]*b*h+c[i,2]*((d*h)+1)+h*exp(logomega[i+1,2])
    logyt[i+1,1:2]~dmnorm(log(c[i+1,1:2]), prec)
    }
    a<- -(beta)/VN
    b<-(beta)/VN
    cc<-beta/VF
    d<- -(beta+Qprime)/VF+kl
    g1<-Gprime/VN
    c[1,1]<-1
    c[1,2]<-1
    logyt[1,1:2]~dmnorm(log(c[1,1:2]), prec)
    prec<-I
    vari~dwish(I,2)
    Qprime~dunif(11,17)
    Gprime~dunif(281,482)
    kl~dunif(0,0.8)
    beta~dunif(0,10)
    }",file="modeljags2zone.jag")
jags<-jags.model("modeljags2zone.jag",data=list("logyt"=logyt, N=99,"I"=diag(2), "VN"=1.1, "VF"=VF,"h"=h))
update(jags, 1000)
mcmcjags.2zone<-jags.samples(jags,
                             c('Qprime','Gprime','beta',"kl",'c'),
                             1000)
mcmcsamplesjags.2zone<-coda.samples(jags,
                                    c('Qprime','Gprime','beta',"kl",'c'),
                                    1000)

quartz()
plot(mcmcsamplesjags.2zone[[1]][,1])
summary(mcmcsamplesjags.2zone)
m1.mcmc.2zone<-(as.mcmc(mcmcsamplesjags.2zone))
m1.mat.2zone<-as.matrix(mcmcsamplesjags.2zone)
m1.dat.2zone<-as.data.frame(m1.mat.2zone)
aomega.post.2zone<-m1.dat.2zone$aomega
bomega.post.2zone<-m1.dat.2zone$bomega
qprime.post.2zone<-m1.dat.2zone$Qprime
gprime.post.2zone<-m1.dat.2zone$Gprime
beta.post.2zone<-m1.dat.2zone$beta
kl.post.2zone<-m1.dat.2zone$kl

c1.hat.2zone <- apply(mcmcsamplesjags.2zone[[1]][,4:103], 2, mean)
c2.hat.2zone <- apply(mcmcsamplesjags.2zone[[1]][,104:203], 2, mean)

c1.hatlower.2zone<- apply(mcmcsamplesjags.2zone[[1]][,4:103], 2, function(x) quantile(x,0.025))
c1.hatupper.2zone<-apply(mcmcsamplesjags.2zone[[1]][,4:103], 2, function(x) quantile(x,0.975))
c2.hatlower.2zone<- apply(mcmcsamplesjags.2zone[[1]][,104:203], 2, function(x) quantile(x,0.025))
c2.hatupper.2zone<-apply(mcmcsamplesjags.2zone[[1]][,104:203], 2, function(x) quantile(x,0.975))
quartz()
plotCI(c[,1], c1.hat.2zone,li=c1.hatlower.2zone, ui=c1.hatupper.2zone)
abline(c(0,1))
plotCI(c[,2], c2.hat.2zone,li=c2.hatlower.2zone, ui=c2.hatupper.2zone)
abline(c(0,1))

pf.post2zone<-A.post<-g.post<-vector("list")
for(i in 1:150){
  A.post[[i]]<-matrix(c(-(beta.post.2zone[i])/VN, (beta.post.2zone[i])/VN, beta.post.2zone[i]/VF, -(beta.post.2zone[i]+qprime.post.2zone[i])/VF+kl.post.2zone[i]),nrow=2, ncol=2,
                      byrow=TRUE)
  g.post[[i]]<-matrix(c(gprime.post.2zone[i]/VN,0),nrow=1, ncol=2)
  pf.post2zone[[i]]<-pf2zone(200,100,logyt,A.post[[i]],g.post[[i]],aomega,bomega,Sigma)
}
logyhat.2zone1<-matrix(0,nrow=100,150)
logyhat.2zone2<-matrix(0,nrow=100,150)
for(i in 1:150){
  for(j in 1:100){
    logyhat.2zone1[,i]<-(pf.post2zone[[i]][1:100])
    logyhat.2zone2[,i]<-(pf.post2zone[[i]][101:200])
    
  }
}
logyhat.2zone1.mean<-apply(logyhat.2zone1,1,mean)
logyhat.2zone2.mean<-apply(logyhat.2zone2,1,mean)

logyt.hatlower.2zone1<- apply(logyhat.2zone1, 1, function(x) quantile(x,0.025))
logyt.hatupper.2zone1<-apply(logyhat.2zone1, 1, function(x) quantile(x,0.975))
logyt.hatlower.2zone2<- apply(logyhat.2zone2, 1, function(x) quantile(x,0.025))
logyt.hatupper.2zone2<-apply(logyhat.2zone2, 1, function(x) quantile(x,0.975))
quartz()
plotCI(truepf.2zone[1:100],logyhat.2zone1.mean[1:100], ui=logyt.hatupper.2zone1[1:100], li=logyt.hatlower.2zone1[1:100])
abline(c(0,1))
quartz()
plotCI(truepf.2zone[101:200],logyhat.2zone2.mean, ui=logyt.hatupper.2zone2, li=logyt.hatlower.2zone2)
abline(c(0,1))
quartz()
points(c[1:100,1],col="blue")
points(c1.hat.2zone, col="yellow")
points(exp(logyt[1:100,1]),col="green")
legend("bottomright",col=c("red","black","blue","yellow","green"),text.font=4,legend=
         c("Particle filter trueN","Particle filter estN","CtN","chatN","ytN"),bg="white", lty=1,cex=0.8)

quartz()
plot(c[,1],col="blue")
points(c1.hat.2zone, col="red")
points(exp(logyt[,1]),col="green")
legend("bottomright",col=c("blue","red","green"),text.font=4,legend=
         c("CtN","chatN","ytN"),bg="white", lty=1,cex=0.8)
#plot(truepf.2zone[101:200], col="red")
#points((logyhat.2zone2.mean))
quartz()
plot(c[,2],col="blue")
points(c2.hat.2zone, col="red")
points(exp(logyt[,2]),col="green")
legend("bottomright",col=c("blue","red","green"),text.font=4,legend=
         c("CtF","chatF","ytF"),bg="white", lty=1,cex=0.8)

theta2julia<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmc2.csv")
apply(theta2julia, 2, function(x){quantile(x)})
xjulia21<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc21.csv")
xjulia22<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc22.csv")
pred21.julia<-apply(xjulia21, 2, mean)
pred22.julia<-apply(xjulia22, 2, mean)

quartz()
plot(truepf.2zone[1:100], col="red")
points((pred21.julia))
points(c[,1],col="blue")
points(exp(logyt[,1]),col="green")
legend("bottomright",col=c("red","black","blue","green"),text.font=4,legend=
         c("Particle filter trueN","Particle filter estN","CtN","ytN"),bg="white", lty=1,cex=0.8)

quartz()
plot(truepf.2zone[101:200], col="red")
points((pred22.julia))
points(c[,2],col="blue")
points(exp(logyt[,2]),col="green")
legend("bottomright",col=c("red","black","blue","green"),text.font=4,legend=
         c("Particle filter trueF","Particle filter estF","CtF","ytF"),bg="white", lty=1,cex=0.8)

#Turbulent eddy-diffusion model
Dt=1
nloc<-5
ntime<-50
x<-seq(1,nloc,length.out=50)
y<-seq(1,nloc,length.out=50)
matdist=NULL
#matdist <- as.matrix(dist(c(x,y),upper=TRUE, diag=TRUE))
for (i in 1:nloc){
  matdist[i]=sqrt(y[i]^2+x[i]^2)
}
erf<-function(x){
  (2/sqrt(pi))*exp(-x^2)
}
csp<-matrix(0,nloc,ntime)
for(i in 1:ntime){
  for(j in 1:nloc){
    set.seed(0000)
    wt<-rgamma(50,1,2)
    csp[j,i]<-Gprime/(2*pi*Dt*(matdist[j]))*(1-
                                               (integrate(erf,lower=0,upper=(matdist[j])/sqrt(4*Dt*i)))$value)+wt[i]
  } 
}

quartz()
plot(csp[1,])
points(csp[2,])
points(csp[3,])
points(csp[4,])
points(csp[5,])


csp2<-matrix(0,nloc,ntime)
csp2[,1]<-csp[,1]
for(i in 2:ntime){
  for(j in 1:nloc){
    csp2[,i]<-csp2[,(i-1)]+h*(Gprime/(4*(pi*Dt*i*h)^1.5))*
      (exp(-(matdist[j])^2/(4*Dt*i*h)))
  }
}
csp2
plot(csp[1,], csp2[1,])
abline(c(0,1))
h=0.5


#quartz()
plot(csp2[1,], ylim=c(0,max(csp2[1,])))
points(csp2[2,])
points(csp2[3,])
points(csp2[4,])
points(csp2[5,])

phi<-1
sigma<-1
set.seed(123)
eta<-rnorm(50,0,0.5)
vt<-NULL
vt[1]<-0
for(i in 2:ntime){
  vt[i]<-vt[i-1]+eta[i]
}
vt
#library(geoR)
Sigma.exp <- function(x,y,phi,sigma) {
  Sigma <- matrix(rep(0, length(x)*length(y)), nrow=length(x))
  Sigma<- as.matrix(exp(-sigma*(dist(data.frame(cbind(x,y)),upper=TRUE, diag=TRUE,method="euclidean")*phi)))
  return(Sigma)
}
sigma.sp<-Sigma.exp(x,y,1,1)
dim(sigma.sp)
set.seed(123)
et<-rnorm(5,rep(0,5),c((0.25*diag(5)+sigma.sp[1:5,1:5])[,1]))
et
# et<-grf(5, grid = "irreg", cbind(x,y), nsim = 1, cov.model = "exponential",
#      cov.pars = c(1,1),
#      kappa = 0.5, nugget = 0, lambda = 1, 
#      mean = 0,  RF=TRUE)$data
logyt3<-matrix(0,nrow=5, ncol=50)
for(j in 1:5){
  logyt3[j,]<-log(csp[j,]+vt)
}
for(i in 1:50){
  logyt3[,i]<-log(exp(logyt3[,i])+et)
}
exp(logyt3)
quartz()
plot(csp[1,])
points(exp(logyt3[1,]),col="red")
quartz()
plot(csp[2,], ylim=c(min(csp[2,]),max(exp(logyt3[2,]))))
points(exp(logyt3[2,]),col="red")
quartz()
plot(csp[3,],ylim=c(min(csp[3,]),max(exp(logyt3[3,]))))
points(exp(logyt3[3,]),col="red")
quartz()
plot(csp[4,])
points(exp(logyt3[4,]),col="red")
quartz()
plot(csp[5,])
points(exp(logyt3[5,]),col="red")

distance<-as.matrix(dist(data.frame(cbind(x,y)),upper=TRUE, diag=TRUE,method="euclidean"))

cat("model{
    for (i in 1:N){
    omega[i+1]~dgamma(aomega,bomega)
    for(j in 1:nloc){
    csp[j,i+1]<-csp[j,i]+h*(Gprime/(4*(pi*Dt*i*h)^1.5))*
    (exp(-(matdist[j])^2/(4*Dt*i*h)))+omega[i+1]
    logyt[j,i+1]~dnorm(log(csp[j,i+1]+et[j]+vt[(i+1)]), tausq)
    }
    }
    for(i in 1:N){
    eta[i]~dnorm(0,2)
    vt[(i+1)]<-vt[i]+eta[i]
    }
    for(j in 1:nloc){
    for(m in 1:nloc){
    k[m,j]<-sigma*exp(-phi*dist[m,j])
    kinv[m,j]<-inverse(k[m,j])
    }
    }
    for(j in 1:nloc){
    et[j]~dnorm(muet,kinv[j,1])
    #logyt[j,1]~dnorm(log(csp[j,1]+et[j]+vt[1]), tausq)
    csp[j,1]<-exp(logyt[j,1])
    }
    vt[1]<-0
    Gprime~dunif(281,482)
    Dt~dunif(0,3)
    sigma~dgamma(2,1)
    phi~dunif(0.5,3)
    aomega~dunif(1,5)
    bomega~dunif(0.5,2)
    }",file="modeljagseddy.jag")
jags.eddy<-jags.model("modeljagseddy.jag",data=list("logyt"=logyt3, N=49,nloc=5,"h"=h,"muet"=0,
                                                    "tausq"=1/0.25,"pi"=pi, "dist"=distance, "matdist"=matdist))
update(jags.eddy, 1000,n.burnin=floor(n.iter/2))
mcmcjags.eddy<-jags.samples(jags.eddy,
                            c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                            1000,n.burnin=floor(n.iter/2))
mcmcsamplesjags.eddy<-coda.samples(jags.eddy,
                                   c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                                   1000,n.burnin=floor(n.iter/2))

quartz()
plot(mcmcsamplesjags.eddy[[1]][,1])
summary(mcmcsamplesjags.eddy)
m1.mcmc.eddy<-(as.mcmc(mcmcsamplesjags.eddy))
m1.mat.eddy<-as.matrix(mcmcsamplesjags.eddy)
m1.dat.eddy<-as.data.frame(m1.mat.eddy)
aomega.post.eddy<-m1.dat.eddy$aomega
bomega.post.eddy<-m1.dat.eddy$bomega
Dt.post.eddy<-m1.dat.eddy$Dt
gprime.post.eddy<-m1.dat.eddy$Gprime
sigma.post.eddy<-m1.dat.eddy$sigma
phi.post.eddy<-m1.dat.eddy$phi
quantile(bomega.post.eddy)

csp.hat<-matrix(0, nrow=5, ncol=50)
for(j in 1:5){
  for(i in 1:50){
    csp.hat[j,i] <- mean(mcmcsamplesjags.eddy[[1]][,4+(5*i+1-(6-j))])
  }
}
csp.hat
plot(csp2[5,],csp.hat[5,])
abline(0,1)
quartz()
plot(csp.hat[1,])
points(csp2[1,],col="blue")
points(exp(logyt3[1,]),col="green")
legend("bottomright",col=c("black","blue","green"),text.font=4,legend=
         c("chat1","Ct1","yt1"),bg="white", lty=1,cex=0.8)

plot(csp.hat[2,])
points(csp2[2,],col="blue")
points(exp(logyt3[2,]),col="green")
plot(csp.hat[3,])
points(csp2[3,],col="blue")
points(exp(logyt3[3,]),col="green")
plot(csp.hat[4,])
points(csp2[4,],col="blue")
points(exp(logyt3[4,]),col="green")
quartz()
plot(csp.hat[5,])
points(csp2[5,],col="blue")
points(exp(logyt3[5,]),col="green")
legend("bottomright",col=c("black","blue","green"),text.font=4,legend=
         c("chat5","Ct5","yt5"),bg="white", lty=1,cex=0.8)
