library(nimble)
library(rjags)
library(plotrix)
library(compositions)
library(mvtnorm)
library(expm)
library(DPpackage)
library(ggvis)
n=100
Qprime=13.8
Gprime=351.54
V=1
sigma<-0.1
aomega<-2
bomega<-1
kl=0.1
set.seed(2017)
wt<-rgamma(n,aomega,bomega)
yt<-NULL
c<-NULL
h=0.01
for(i in 1:(n-1)){
  c[1]<-1+wt[1]
  c[i+1]<-(1-h*(Qprime+kl*V)/V)*c[i]+h*Gprime/V+wt[i+1]-wt[i]
}
#exact solution
c2<-NULL
for(i in 1:(n-1)){
  c2[1]<-1+wt[1]
  c2[i+1]<-exp(-h*i*(Qprime+kl*V)/V)*c[1]+1/(-(Qprime+kl*V)/V)*
    (exp(-h*i*(Qprime+kl*V)/V)-1)*Gprime/V+wt[i+1]
}
c2
plot(c2)
plot(c2,c)
abline(0,1)
Qprime/V+kl
set.seed(000)
vt<-rnorm(n,0,sigma)

logyt<-log(c2)+(vt)
plot(c2)
points(exp(logyt),col="green")
hist(exp(logyt))

#jags
cat("model{
    for (i in 1:N){  
    omega[i+1]~dgamma(aomega,bomega)
    c[i+1]<-(1-h*(Qprime+kl*V)/V)*c[i]+h*Gprime/V+(omega[i+1]-omega[i])
    logyt[i+1]~dnorm(log(c[i+1])+muv[z1[i+1]],sigmav[z2[i+1]])
    z1[i+1]~dcat(p1[])
    z2[i+1]~dcat(p2[])
    }
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
    basetau<-0.1
    omega[1]~dgamma(aomega,bomega)
    c[1]<-1
    logyt[1]~dnorm(log(c[1]), prec)
    alpha<-1
    aomega~dunif(0.5,3)
    bomega~dunif(0.5,3)
    sigma~dgamma(2,0.01)
    prec<-1/sigma
    Qprime~dunif(11,17)
    Gprime~dunif(281,482)
    kl~dunif(0,0.8)
    }",file="modeljagsdp.jag")
jags.dp<-jags.model("modeljagsdp.jag",data=list("logyt"=logyt, N=99,"h"=0.01, "V"=1,"m"=100))
update(jags.dp, 1000)
mcmcjags.dp<-jags.samples(jags.dp,
                          c('Qprime','Gprime','kl','sigma','c','aomega','bomega','alpha'),
                          1000)
mcmcsamplesjags.dp<-coda.samples(jags.dp,
                                 c('Qprime','Gprime','kl','sigma','c','aomega','bomega','alpha'),
                                 1000)
quartz()
plot(mcmcsamplesjags.dp)
m1.mcmc.dp<-(as.mcmc(mcmcsamplesjags.dp))
m1.mat.dp<-as.matrix(mcmcsamplesjags.dp)
m1.dat.dp<-as.data.frame(m1.mat.dp)
aomega.post1zone.dp<-m1.dat.dp$aomega
bomega.post1zone.dp<-m1.dat.dp$bomega
qprime.post1zone.dp<-m1.dat.dp$Qprime
gprime.post1zone.dp<-m1.dat.dp$Gprime
kl.post1zone.dp<-m1.dat.dp$kl
quantile(qprime.post1zone.dp,c(0.025,0.5,0.975))
quantile(gprime.post1zone.dp,c(0.025,0.5,0.975))
quantile(kl.post1zone.dp,c(0.025,0.5,0.975))
quantile(aomega.post1zone.dp,c(0.025,0.5,0.975))
quantile(bomega.post1zone.dp,c(0.025,0.5,0.975))
c.hat.dp <- apply(mcmcjags.dp$c, 1, mean)
c.hatlower.dp<- apply(mcmcjags.dp$c, 1, function(x) quantile(x,0.025))
c.hatupper.dp<-apply(mcmcjags.dp$c, 1, function(x) quantile(x,0.975))
quartz()
plotCI(c, c.hat.dp,ui=c.hatupper.dp, li=c.hatlower.dp)
abline(c(0,1))

quartz()
plot(c,col="blue")
points(c.hat.dp, col="red")
points(exp(logyt),col="green")
legend("bottomrigh", legend=c("True","Predicted","Observed"),col=c("blue","red","green"),
       pch=15,cex=1)
#interactive plot
d1<-data.frame(cbind(t=seq(1,100,1),c,c.hat.dp,y=exp(logyt)))
quartz()
d1%>%ggvis(~t,~c)%>%
  layer_points(fill:=input_select(c("red","blue","green")))%>%
  layer_points(x=~t, y=~c.hat.dp,fill:=input_select(c("red","blue","green")))%>%
  layer_points(x=~t, y=~y,fill:=input_select(c("red","blue","green")))%>%
  add_axis("x",title="Time")%>%
  add_axis("y",title="Concentrations")

#DP all
#jags
cat("model{
    for (i in 1:N){ 
    omega[i+1]~dnorm(aomega[z1[i+1]],bomega[z2[i+1]])
    c[i+1]<-(1-h*(Qprime+kl*V)/V)*c[i]+h*Gprime/V+(omega[i+1])
    logyt[i+1]~dnorm(log(c[i+1])+muv[z3[i+1]],sigmav[z4[i+1]])
    z1[i+1]~dcat(p1[])
    z2[i+1]~dcat(p2[])
    z3[i+1]~dcat(p3[])
    z4[i+1]~dcat(p4[])    
    }
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
    r3[l]~dbeta(1,alpha1)
    r4[l]~dbeta(1,alpha1)
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
    aomega[l]~dnorm(basemu,basetau)
    bomega[l]~dgamma(2,1)
    #aomega[l]~dunif(1,5)
    #bomega[l]~dunif(0.5,2)
    }
    basemu<-0
    basetau<-1
    omega[1]~dnorm(aomega[z1[2]],bomega[z2[2]])
    c[1]<-1
    logyt[1]~dnorm(log(c[1]), sigmav[1])
    alpha<-1
    alpha1<-1
    Qprime~dunif(11,17)
    Gprime~dunif(281,482)
    kl~dunif(0,0.8)
    }",file="modeljagsdp2.jag")

jags.dp2<-jags.model("modeljagsdp2.jag",data=list("logyt"=logyt, N=99, "h"=0.01,"V"=1,"m"=100))
update(jags.dp2, 1000)
mcmcjags.dp2<-jags.samples(jags.dp2,
                           c('Qprime','Gprime','kl','sigmav','c','aomega','bomega'),
                           1000)
mcmcsamplesjags.dp2<-coda.samples(jags.dp2,
                                  c('Qprime','Gprime','kl','sigmav','c','aomega','bomega'),
                                  1000)
m1.mcmc.dp2<-(as.mcmc(mcmcsamplesjags.dp2))
m1.mat.dp2<-as.matrix(mcmcsamplesjags.dp2)
m1.dat.dp2<-as.data.frame(m1.mat.dp2)
qprime.post1zone.dp2<-m1.dat.dp2$Qprime
gprime.post1zone.dp2<-m1.dat.dp2$Gprime
kl.post1zone.dp2<-m1.dat.dp2$kl
quantile(qprime.post1zone.dp2,c(0.025,0.5,0.975))
quantile(gprime.post1zone.dp2,c(0.025,0.5,0.975))
quantile(kl.post1zone.dp2,c(0.025,0.5,0.975))

c.hat.dp2 <- apply(mcmcjags.dp2$c, 1, mean)
c.hatlowerdp2<- apply(mcmcjags.dp2$c, 1, function(x) quantile(x,0.025))
c.hatupperdp2<-apply(mcmcjags.dp2$c, 1, function(x) quantile(x,0.975))
quartz()
plotCI(c, c.hat.dp2[1:100],ui=c.hatupperdp2[1:100], li=c.hatlowerdp2[1:100])
abline(0,1)
quartz()
plot((c),col="blue")
points(c.hat.dp2,col="red")
points(exp(logyt), col="green")
legend("bottomright",legend=c("True","Predicted","Observed"),cex=1, pch=15,
       col=c("blue","red","green"))
