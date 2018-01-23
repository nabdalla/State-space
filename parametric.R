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
#values of c very huge or very small and hence yt so taking log is problem
for(i in 1:(n-1)){
  c[1]<-1+wt[1]
  c[i+1]<-(1-h*(Qprime+kl*V)/V)*c[i]+h*Gprime/V+wt[i+1]-wt[i]
}
plot(c)
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
    logyt[i+1]~dnorm(log(c[i+1]), prec)
    }
    omega[1]~dgamma(aomega,bomega)
    c[1]<-1
    logyt[1]~dnorm(log(c[1])+omega[1], prec)
    sigma<-1/prec
    aomega~dunif(0.5,3)
    bomega~dunif(0.5,3)
    prec~dgamma(2,0.01)
    Qprime~dunif(11,17)
    Gprime~dunif(281,482)
    kl~dunif(0,0.8)
    for(i in (N+1):(N+9)){
    omega[i+1]~dgamma(aomega,bomega)
    c[i+1]<-(1-h*(Qprime+kl*V)/V)*c[i]+h*Gprime/V+(omega[i+1]-omega[i])
    pred[i-N]~dnorm(log(c[i]),prec)
    }
    }",file="modeljags.jag")
jags<-jags.model("modeljags.jag",data=list("logyt"=logyt[1:100], N=99, "V"=V,"h"=h))
update(jags, 1000)
mcmcjags<-jags.samples(jags,
                       c('Qprime','Gprime',"kl",'sigma','c','aomega','bomega','pred'),
                       1000)
mcmcsamplesjags<-coda.samples(jags,
                              c('Qprime','Gprime','kl','sigma','c','aomega','bomega','pred'),
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
sigma.post1zone<-sqrt(m1.dat$sigma)
pred.post1zone<-m1.dat$pred
quantile(sigma.post1zone,c(0.025,0.25,0.5,0.75,0.975))
c.hat <- apply(mcmcjags$c, 1, mean)
c.hatlower<- apply(mcmcjags$c, 1, function(x) quantile(x,0.025))
c.hatupper<-apply(mcmcjags$c, 1, function(x) quantile(x,0.975))
pred.hat<-apply(mcmcjags$pred,1,mean)
plot(logyt[101:109])
points(pred.hat,col="red")
abline(c(0,1))
quartz()
plotCI(c2, c.hat,ui=c.hatupper, li=c.hatlower)
abline(c(0,1))
quartz()
plot(exp(logyt),col="green")
points(c2,col="blue")
points(c.hat,col="red")
mse1sim<-sum((exp(logyt[101:109])-exp(pred.hat))^2)/length(pred.hat)
#Pf import from julia output
thetajulia<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmc.csv")
apply(thetajulia, 2, function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})
plot(thetajulia[,3],type="l")
acf(thetajulia[,2], lag.max=1000)
xjulia<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc.csv")
pred.julia<-apply(xjulia, 2, mean)
quartz()
plot(exp(logyt),col="green")
points(c2,col="blue")
points(pred.julia,col="red")
#nonlinear regression
cat("model{
    for (i in 1:N){ 
    c[i+1]<-exp(-i*h*(Qprime+kl*V)/V)*c[1]+1/(-(Qprime+kl*V)/V)*
    (exp(-i*h*(Qprime+kl*V)/V)-1)*Gprime/V
    logyt[i+1]~dnorm(log(c[i+1]), prec)
    }
    c[1]<-1
    logyt[1]~dnorm(log(c[1]), prec)
    sigma<-1/prec
    prec~dgamma(2,0.01)
    Qprime~dunif(11,17)
    Gprime~dunif(281,482)
    kl~dunif(0,0.8)
    }",file="modeljagsnl.jag")
jagsnl<-jags.model("modeljagsnl.jag",data=list("logyt"=logyt[1:100], N=99, "V"=V,"h"=0.01))
update(jags, 1000)
mcmcjagsnl<-jags.samples(jagsnl,
                         c('Qprime','Gprime',"kl",'sigma','c'),
                         1000)
mcmcsamplesjagsnl<-coda.samples(jagsnl,
                                c('Qprime','Gprime','kl','sigma','c'),
                                1000)

summary(mcmcsamplesjagsnl)
m1.mcmcnl<-(as.mcmc(mcmcsamplesjagsnl))
m1.matnl<-as.matrix(mcmcsamplesjagsnl)
m1.datnl<-as.data.frame(m1.matnl)
c.hatnl <- apply(mcmcjagsnl$c, 1, mean)

plot(exp(logyt),col="green")
points(c2,col="blue")
points(c.hatnl,col="red")

#all simulation plots
quartz()

layout(matrix(c(1,1,1,1,1,0,2,2,2,2,2,0,0,0,3,3,3,3,3,
                0,0,0), 2, 11, byrow = TRUE))
plot(exp(logyt),col="green", main="State-by-state update", ylab="Concentrations", xlab="Time")
points(c2,col="blue")
points(c.hat,col="red")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)
plot(exp(logyt), col="green",main="KF",ylab="Concentrations", xlab="Time")
points(predkf,col="red")
points(c2,col="blue")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)

plot(exp(logyt),col="green", main="PMMH",ylab="Concentrations", xlab="Time")
points(c2,col="blue")
points(pred.julia,col="red")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)
#data one compartment
data1<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/Modeling_Master_v5.csv")
dim(data1)
data1zone=matrix(0,nrow=204, ncol=(81*6))
for(i in 1:(81*6)){
  data1zone[,i]=data1[,4*i]
}
data1zone[,204]
write.csv(data1zone,"/Users/n_a_abdallah/Desktop/spatial/Project2/acetone1zone.csv")
#true: V=11.9, Q=0.04-0.07, 0.23-0.27 and 0.47-0.77 m^3/min, 
#G=39.6, 79.1 and 118.7 for l, med, h
#note that each was measured at 6 locations
quartz()
plot(data1zone[,2])
points(data1zone[,30])
points(data1zone[,51*4])
points(data1zone[,400])
points(data1zone[,55])
points(data1zone[,6])

cat("model{
    for (i in 1:N){ 
    omega[i+1]~dgamma(aomega,bomega)
    omega1[i+1]~dgamma(aomega,bomega)
    c[i+1]<-(1-h*(Qprime+kl*V)/V)*c[i]+h*Gprime/V+omega[i+1]-omega1[i+1]
    logyt[i+1]~dnorm(log(abs(c[i+1])), prec)
    }
    omega[1]~dgamma(aomega,bomega)
    c[1]<-0+omega[1]
    logyt[1]~dnorm(log(c[1]), prec)
    sigma<-1/prec
    aomega~dunif(0.5,3)
    bomega~dunif(0.5,3)
    prec~dgamma(1,0.01)
    Qprime~dunif(0,1)
    Gprime~dunif(30,150)
    kl~dunif(0,0.8)
    }",file="modeljags.jag")
jagslow<-jags.model("modeljags.jag",data=list("logyt"=log(data1zone[1:201,2]), N=200, "V"=11.9,"h"=0.01))
update(jagslow, 1000)
mcmcjagslow<-jags.samples(jagslow,
                          c('Qprime','Gprime',"kl",'sigma','c','aomega','bomega'),
                          1000)
mcmcsamplesjagslow<-coda.samples(jagslow,
                                 c('Qprime','Gprime','kl','sigma','c','aomega','bomega'),
                                 1000)
#low ventilation Q=0.04-0.07 G=43.18
quartz()
plot(mcmcsamplesjagslow)
summary(mcmcsamplesjagslow)
m1.mcmclow<-(as.mcmc(mcmcsamplesjagslow))
m1.matlow<-as.matrix(mcmcsamplesjagslow)
m1.datlow<-as.data.frame(m1.matlow)
aomega.post1zonelow<-m1.datlow$aomega
bomega.post1zonelow<-m1.datlow$bomega
qprime.post1zonelow<-m1.datlow$Qprime
gprime.post1zonelow<-m1.datlow$Gprime
kl.post1zonelow<-m1.datlow$kl
quantile((gprime.post1zonelow),c(0.025,0.25,0.5,0.75,0.975))
quantile((qprime.post1zonelow),c(0.025,0.25,0.5,0.75,0.975))
quantile((kl.post1zonelow),c(0.025,0.25,0.5,0.75,0.975))
c.hatlow <- apply(mcmcjagslow$c, 1, mean)

#medium Q Q=0.23-0.27 G=43.18
jagsmed<-jags.model("modeljags.jag",data=list("logyt"=log(data1zone[,30]+0.05), N=80, "V"=11.9,"h"=0.01))
update(jagsmed, 1000)
mcmcjagsmed<-jags.samples(jagsmed,
                          c('Qprime','Gprime',"kl",'sigma','c','aomega','bomega'),
                          1000)
mcmcsamplesjagsmed<-coda.samples(jagsmed,
                                 c('Qprime','Gprime','kl','sigma','c','aomega','bomega'),
                                 1000)
summary(mcmcsamplesjagsmed)
m1.mcmcmed<-(as.mcmc(mcmcsamplesjagsmed))
m1.matmed<-as.matrix(mcmcsamplesjagsmed)
m1.datmed<-as.data.frame(m1.matmed)
gprime.post1zonemed<-m1.datmed$Gprime
qprime.post1zonemed<-m1.datmed$Qprime
quantile((gprime.post1zonemed),c(0.025,0.25,0.5,0.75,0.975))
quantile((qprime.post1zonemed),c(0.025,0.25,0.5,0.75,0.975))
c.hatmed <- apply(mcmcjagsmed$c, 1, mean)

#high ventilation Q=0.47-0.77 low G=39.55
jagsh<-jags.model("modeljags.jag",data=list("logyt"=log(data1zone[,51*4]+0.05), N=40, "V"=11.9,"h"=0.01))
update(jagsh, 1000)
mcmcjagsh<-jags.samples(jagsh,
                        c('Qprime','Gprime',"kl",'sigma','c','aomega','bomega'),
                        1000)
mcmcsamplesjagsh<-coda.samples(jagsh,
                               c('Qprime','Gprime','kl','sigma','c','aomega','bomega'),
                               1000)
summary(mcmcsamplesjagsh)
m1.mcmch<-(as.mcmc(mcmcsamplesjagsh))
m1.math<-as.matrix(mcmcsamplesjagsh)
m1.dath<-as.data.frame(m1.math)
gprime.post1zoneh<-m1.dath$Gprime
qprime.post1zoneh<-m1.dath$Qprime
quantile((gprime.post1zoneh),c(0.025,0.25,0.5,0.75,0.975))
quantile((qprime.post1zoneh),c(0.025,0.25,0.5,0.75,0.975))
c.hath <- apply(mcmcjagsh$c, 1, mean)
quartz()
par(mfrow=c(1,3))
plot(c.hatlow, main="low ventilation")
points(data1zone[1:201,2],col="red")
legend("bottomright", col=c("black","red"),legend=c("chat low","measured"))
plot(c.hatmed, main="medium ventilation")
points(data1zone[1:81,30],col="red")
legend("bottomright", col=c("black","red"),legend=c("chat medium","measured"))
plot(c.hath,main="High ventilation")
points(data1zone[,51*4],col="red")
legend("bottomright", col=c("black","red"),legend=c("chat high","measured"))
#Pf data one zone
#Pf import from julia output
#low
thetajuliadl<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmcl.csv")
apply(thetajuliadl, 2, function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})
xjuliadl<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmcl.csv")
pred.juliadl<-apply(xjuliadl, 2, mean)

#medium
thetajuliadm<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmcm.csv")
apply(thetajuliadm, 2, function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})
xjuliadm<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmcm.csv")
pred.juliadm<-apply(xjuliadm, 2, mean)

#high
thetajuliadh<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmch.csv")
apply(thetajuliadh, 2, function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975))})
xjuliadh<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmch.csv")
pred.juliadh<-apply(xjuliadh, 2, mean)
#plot all data
quartz()
par(mfrow=c(1,3))
plot(data1zone[1:201,2], col="green", main="Low", ylab="Toluene Concentrations",
     xlab="Time")
points(pred.juliadl,col="red")
points(c.hatlow, col="blue")
legend("bottomright", col=c("green", "blue", "red"), legend=c("Measurements",
                                                              "State-by-state update", "PMMH"),pch=15,cex=1,box.lwd = 0)
plot(pred.juliadm,col="red", main="Medium", ylab="Toluene Concentrations",
     xlab="Time")
points(data1zone[1:81,30], col="green")
points(c.hatmed, col="blue")
legend("bottomright", col=c("green", "blue", "red"), legend=c("Measurements",
                                                              "State-by-state update", "PMMH"),pch=15,cex=1,box.lwd = 0)
plot(pred.juliadh,col="red", main="High", ylab="Acetone Concentrations",
     xlab="Time")
points(data1zone[1:41,203], col="green")
points(c.hath, col="blue")
legend("bottomright", col=c("green", "blue", "red"), legend=c("Measurements",
                                                              "State-by-state update", "PMMH"),pch=15,cex=1,box.lwd = 0)
#two compartment model simulation
n=100
Qprime=13.8
Gprime=351.54
VF=1
VN=1
Sigma<-diag(0.1,2)
kl=0.1
beta=5
V=diag(0.1,2)
r=3
set.seed(2017)
wt<-rlnorm.rplus(n,c(0,0),diag(2))
#wt<-rmvnorm(n,c(0,0),V)
h=0.01
A<-matrix(c(-(beta)/VN, (beta)/VN, beta/VF, -(beta+Qprime)/VF+kl),nrow=2, ncol=2,
          byrow=TRUE)
eigen(A)
g<-matrix(c(Gprime/VN,0),nrow=1, ncol=2)

yt<-NULL
c<-matrix(c(rep(0,2*n)), nrow=n, ncol=2)
#values of c very huge or very small and hence yt so taking log is problem
for(i in 1:(n-1)){
  c[1,]<-c(0,0.5)+wt[1,]
  c[i+1,]<-(c[i,])%*%(h*A+diag(1,2))+(g*h)+h*wt[i+1,]
}
c
quartz()
plot(c[,1],ylab="mg/m3", xlab="time")
points(c[,2])

#exact 
c2<-matrix(c(rep(0,2*n)), nrow=n, ncol=2)
#imp

for(i in 1:(n-1)){
  c2[1,]<-c(0,0.5)+wt[1,]
  c2[i+1,]<- expm(h*i*A)%*%c2[1,]+solve(A)%*%(expm(h*i*A)-diag(2))%*%t(g)+wt[i+1,]
}
plot(c2[,1])
points(c2[,2])
plot(c[,1]/VN,c2[,1])
abline(0,1)
plot(c[,2]/VF,c2[,2])
abline(0,1)

c<-c2

set.seed(123)
vt<-rmvnorm(n,c(0,0),Sigma)
logyt2<-log(abs(c))+(vt)
exp(logyt2)

#jags
I<-diag(100,2)

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
    prec~dwish(I,2)
    vari~dwish(I,2)
    Qprime~dunif(11,17)
    Gprime~dunif(281,482)
    kl~dunif(0,0.8)
    beta~dunif(0,10)
    }",file="modeljags2zone.jag")
jags<-jags.model("modeljags2zone.jag",data=list("logyt"=logyt2, N=99,"I"=diag(2), "VN"=1, "VF"=1,"h"=h))
update(jags, 1000)
mcmcjags.2zone<-jags.samples(jags,
                             c('Qprime','Gprime','beta',"kl",'c',"prec","vari"),
                             1000)
mcmcsamplesjags.2zone<-coda.samples(jags,
                                    c('Qprime','Gprime','beta',"kl",'c',"prec","vari"),
                                    1000)

quartz()
plot(mcmcsamplesjags.2zone[[1]][,1])
summary(mcmcsamplesjags.2zone[[1]])
summary(1/mcmcsamplesjags.2zone[[1]][,205:212])
m1.mcmc.2zone<-(as.mcmc(mcmcsamplesjags.2zone))
m1.mat.2zone<-as.matrix(mcmcsamplesjags.2zone)
m1.dat.2zone<-as.data.frame(m1.mat.2zone)
qprime.post.2zone<-m1.dat.2zone$Qprime
gprime.post.2zone<-m1.dat.2zone$Gprime
beta.post.2zone<-m1.dat.2zone$beta
kl.post.2zone<-m1.dat.2zone$kl
quantile(kl.post.2zone,c(0.025,0.25,0.5,0.75,0.975))
quantile(prec.post.2zone,c(0.025,0.25,0.5,0.75,0.975))
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

quartz()
par(mfrow=c(1,2))
plot(c[,1],col="red")
points(c1.hat.2zone, col="blue")
points(exp(logyt2[,1]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("CtN","chatN","ytN"),bg="white", lty=1,cex=0.8)
plot(c[,2],col="red")
points(c2.hat.2zone, col="blue")
points(exp(logyt2[,2]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("CtF","chatF","ytF"),bg="white", lty=1,cex=0.8)

theta2julia<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmc2.csv")
apply(theta2julia, 2, function(x){quantile(x)})
xjulia21<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc21.csv")
xjulia22<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc22.csv")
pred21.julia<-apply(xjulia21, 2, mean)
pred22.julia<-apply(xjulia22, 2, mean)

quartz()
par(mfrow=c(1,2))
plot(c[,1],col="blue", main="State-by-state Gibbs Near", ylab="Concentrations", xlab="Time")
points(pred21.julia, col="red")
points(exp(logyt[,1]),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements N","True N",
                                                            "Estimated N"),pch=15,cex=0.75,box.lwd = 0)
plot(c[,2],col="blue", main="State-by-state Gibbs Far", ylab="Concentrations", xlab="Time")
points(pred22.julia, col="red")
points(exp(logyt[,2]),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements F","True F",
                                                            "Estimated F"),pch=15,cex=0.75,box.lwd = 0)

#all sim plots
quartz()
par(mfrow=c(3,2))
plot(c[,1],col="blue", main="State-by-state update Near", ylab="Concentrations", xlab="Time")
points(c1.hat.2zone, col="red")
points(exp(logyt2[,1]),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements N","True N",
                                                            "Estimated N"),pch=15,cex=0.75,box.lwd = 0)
plot(c[,2],col="blue", main="State-by-state update Far", ylab="Concentrations", xlab="Time")
points(c2.hat.2zone, col="red")
points(exp(logyt2[,2]),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements F","True F",
                                                            "Estimated F"),pch=15,cex=0.75,box.lwd = 0)
plot(exp(logyt2[,1]), col="green", main="KF Near", ylab="Concentrations", xlab="Time")
points(predkf1,col="red")
points(c[,1],col="blue")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements N","True N",
                                                            "Estimated N"),pch=15,cex=0.75,box.lwd = 0)
plot(exp(logyt2[,2]), col="green", main="KF Far", ylab="Concentrations", xlab="Time")
points(predkf2,col="red")
points(c[,2],col="blue")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements F","True F",
                                                            "Estimated F"),pch=15,cex=0.75,box.lwd = 0)
plot(c[,1],col="blue", main="PMMH Near", ylab="Concentrations", xlab="Time")
points(pred21.julia, col="red")
points(exp(logyt2[,1]),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements N","True N",
                                                            "Estimated N"),pch=15,cex=0.75,box.lwd = 0)
plot(c[2:100,2],col="blue", main="PMMH Far", ylab="Concentrations", xlab="Time")
points(pred22.julia[2:100], col="red")
points(exp(logyt2[2:100,2]),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements F","True F",
                                                            "Estimated F"),pch=15,cex=0.75,box.lwd = 0)

#data
data2zone=matrix(0,nrow=204, ncol=(81*4*6))
data2zonen=matrix(0,nrow=204, ncol=81)
data2zonef1=matrix(0,nrow=204, ncol=81)
data2zonef2=matrix(0,nrow=204, ncol=81)
data2zonef3=matrix(0,nrow=204, ncol=81)
dim(data1)
for(i in 1:(81*6)){
  data2zone[,i]=data1[,(81*4*6)+(4*i)]
}
for(i in 1:81){
  data2zonen[,i]=data2zone[,1+(6*(i-1))]
  data2zonef1[,i]=data2zone[,4+(6*(i-1))]
  data2zonef2[,i]=data2zone[,5+(6*(i-1))]
  data2zonef3[,i]=data2zone[,6+(6*(i-1))]
}
write.csv(data2zonen,"/Users/n_a_abdallah/Desktop/spatial/Project2/data2zonen.csv")
write.csv(data2zonef1,"/Users/n_a_abdallah/Desktop/spatial/Project2/data2zonef1.csv")
write.csv(data2zonef2,"/Users/n_a_abdallah/Desktop/spatial/Project2/data2zonef2.csv")
write.csv(data2zonef3,"/Users/n_a_abdallah/Desktop/spatial/Project2/data2zonef3.csv")

#true values Q=0.059, 0.258 and 0.595 Vn=0.105 Vf=11.79 beta=0.24-1.24
quartz()
plot(data2zonef1[,10])
points(data2zonef1[,81])
cat("model{
    for (i in 1:N){
    logomega[i+1,1:2]~dmnorm(c(0,0),vari)
    c[i+1,1]<-c[i,1]*((a*h)+1)+c[i,2]*b*h+g1*h+h*exp(logomega[i+1,1])
    c[i+1,2]<-c[i,1]*cc*h+c[i,2]*((d)+1)+exp(logomega[i+1,2])
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
    prec~dwish(I,10)
    vari~dwish(I,10)
    Qprime~dunif(0,1)
    Gprime~dunif(30,150)
    kl~dunif(0,0.8)
    beta~dunif(0,5)
    }",file="modeljags2zone.jag")
#med ventilation Q=0.23-0.27 med G=86.36
jagsm<-jags.model("modeljags2zone.jag",data=list("logyt"=cbind(log(data2zonen[1:75,13]+0.05),log(data2zonef1[1:75,13]+0.05)),
                                                 N=74,"I"=diag(2), "VN"=0.1, "VF"=11.9,"h"=0.01),
                  n.adapt =1000)
update(jagsm, 2000)
mcmcjags.2zonem<-jags.samples(jagsm,
                              c('Qprime','Gprime','beta',"kl",'c',"prec","vari"),
                              1000)
mcmcsamplesjags.2zonem<-coda.samples(jagsm,
                                     c('Qprime','Gprime','beta',"kl",'c',"prec","vari"),
                                     1000)

quartz()
plot(mcmcsamplesjags.2zonem[[1]][,1])
summary(mcmcsamplesjags.2zonem[[1]])
m1.mcmc.2zonem<-(as.mcmc(mcmcsamplesjags.2zonem))
m1.mat.2zonem<-as.matrix(mcmcsamplesjags.2zonem)
m1.dat.2zonem<-as.data.frame(m1.mat.2zonem)
qprime.post.2zonem<-m1.dat.2zonem$Qprime
gprime.post.2zonem<-m1.dat.2zonem$Gprime
beta.post.2zonem<-m1.dat.2zonem$beta
kl.post.2zonem<-m1.dat.2zonem$kl
quantile(qprime.post.2zonem,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post.2zonem,c(0.025,0.25,0.5,0.75,0.975))
quantile(beta.post.2zonem,c(0.025,0.25,0.5,0.75,0.975))
c1.hat.2zonem <- apply(mcmcsamplesjags.2zonem[[1]][,4:(78)], 2, mean)
c2.hat.2zonem <- apply(mcmcsamplesjags.2zonem[[1]][,(79):(78+75)], 2, mean)
quartz()
plot(c2.hat.2zonem,data2zonef1[1:75,13])
abline(0,1)

quartz()
plot(data2zonen[,3])
points(data2zonef1[,3])
#high
jagsh<-jags.model("modeljags2zone.jag",data=list("logyt"=cbind(log(data2zonen[,81]),log(data2zonef1[,81])),
                                                 N=40,"I"=diag(2), "VN"=0.1, "VF"=11.9,"h"=0.01),
                  n.adapt =1000)
update(jagsh, 2000)
mcmcjags.2zoneh<-jags.samples(jagsh,
                              c('Qprime','Gprime','beta',"kl",'c',"prec","vari"),
                              1000)
mcmcsamplesjags.2zoneh<-coda.samples(jagsh,
                                     c('Qprime','Gprime','beta',"kl",'c',"prec","vari"),
                                     1000)

quartz()
plot(mcmcsamplesjags.2zoneh[[1]][,1])
summary(mcmcsamplesjags.2zoneh[[1]])
m1.mcmc.2zoneh<-(as.mcmc(mcmcsamplesjags.2zoneh))
m1.mat.2zoneh<-as.matrix(mcmcsamplesjags.2zoneh)
m1.dat.2zoneh<-as.data.frame(m1.mat.2zoneh)
qprime.post.2zoneh<-m1.dat.2zoneh$Qprime
gprime.post.2zoneh<-m1.dat.2zoneh$Gprime
beta.post.2zoneh<-m1.dat.2zoneh$beta
kl.post.2zoneh<-m1.dat.2zoneh$kl
quantile(qprime.post.2zoneh,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post.2zoneh,c(0.025,0.25,0.5,0.75,0.975))
quantile(beta.post.2zoneh,c(0.025,0.25,0.5,0.75,0.975))
c1.hat.2zoneh <- apply(mcmcsamplesjags.2zoneh[[1]][,4:(44)], 2, mean)
c2.hat.2zoneh <- apply(mcmcsamplesjags.2zoneh[[1]][,(45):85], 2, mean)
#low Q=0.056
jagsl<-jags.model("modeljags2zone.jag",data=list("logyt"=cbind(log(data2zonen[1:176,1]),log(data2zonef1[1:176,1]+0.05)),
                                                 N=175,"I"=diag(2), "VN"=0.1, "VF"=11.9,"h"=0.01),
                  n.adapt =1000)
update(jagsl, 2000)
mcmcjags.2zonel<-jags.samples(jagsl,
                              c('Qprime','Gprime','beta',"kl",'c',"prec","vari"),
                              1000)
mcmcsamplesjags.2zonel<-coda.samples(jagsl,
                                     c('Qprime','Gprime','beta',"kl",'c',"prec","vari"),
                                     1000)

quartz()
plot(mcmcsamplesjags.2zonel[[1]][,1])
summary(mcmcsamplesjags.2zonel[[1]])
m1.mcmc.2zonel<-(as.mcmc(mcmcsamplesjags.2zonel))
m1.mat.2zonel<-as.matrix(mcmcsamplesjags.2zonel)
m1.dat.2zonel<-as.data.frame(m1.mat.2zonel)
qprime.post.2zonel<-m1.dat.2zonel$Qprime
gprime.post.2zonel<-m1.dat.2zonel$Gprime
beta.post.2zonel<-m1.dat.2zonel$beta
kl.post.2zonel<-m1.dat.2zonel$kl
quantile(qprime.post.2zonel,c(0.025,0.25,0.5,0.75,0.975))
quantile(gprime.post.2zonel,c(0.025,0.25,0.5,0.75,0.975))
quantile(beta.post.2zonel,c(0.025,0.25,0.5,0.75,0.975))
c1.hat.2zonel <- apply(mcmcsamplesjags.2zonel[[1]][,4:(103+76)], 2, mean)
c2.hat.2zonel <- apply(mcmcsamplesjags.2zonel[[1]][,(104+76):355], 2, mean)
quartz()
plot(c2.hat.2zonel,data2zonef1[1:176,1])
abline(0,1)
quartz()
par(mfrow=c(3,2))
plot(c1.hat.2zonel, col="blue",main="Low")
points(data2zonen[1:176,1],col="red")
legend("bottomright",col=c("blue","red"),legend=c("chat near","measured near"),
       bg="white", lty=1,cex=0.8)
plot(c2.hat.2zonel, col="blue")
points(data2zonef1[1:176,1],col="red")
legend("bottomright",col=c("blue","red"),legend=c("chat far","measured far"),
       bg="white", lty=1,cex=0.8)

plot(c1.hat.2zoneh, col="blue",main="Medium")
points(data2zonen[1:75,13],col="red")
legend("bottomright",col=c("blue","red"),legend=c("chat near","measured near"),
       bg="white", lty=1,cex=0.8)
plot(c2.hat.2zoneh, col="blue")
points(data2zonef1[1:75,13],col="red")
legend("bottomright",col=c("blue","red"),legend=c("chat far","measured far"),
       bg="white", lty=1,cex=0.8)

plot(c1.hat.2zonem, col="blue",main="High")
points(data2zonen[1:41,81],col="red")
legend("bottomright",col=c("blue","red"),legend=c("chat near","measured near"),
       bg="white", lty=1,cex=0.8)
plot(c2.hat.2zonem, col="blue")
points(data2zonef1[1:41,81],col="red")
legend("bottomright",col=c("blue","red"),legend=c("chat far","measured far"),
       bg="white", lty=1,cex=0.8)

#PF
theta2julial<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmc2l.csv")
apply(theta2julial, 2, function(x){quantile(x)})
xjulia21l<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc21l.csv")
xjulia22l<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc22l.csv")
pred21.julial<-apply(xjulia21l[2000:5000,], 2, mean)
pred22.julial<-apply(xjulia22l, 2, mean)

theta2juliam<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmc2m.csv")
apply(theta2juliam, 2, function(x){quantile(x)})
xjulia21m<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc21m.csv")
xjulia22m<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc22m.csv")
pred21.juliam<-apply(xjulia21m, 2, mean)
pred22.juliam<-apply(xjulia22m, 2, mean)

theta2juliah<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmc2h.csv")
apply(theta2juliah, 2, function(x){quantile(x)})
xjulia21h<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc21h.csv")
xjulia22h<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc22h.csv")
pred21.juliah<-apply(xjulia21h, 2, mean)
pred22.juliah<-apply(xjulia22h, 2, mean)
quartz()
par(mfrow=c(3,2))
plot(pred21.julial, col="blue",main="Low")
points(data2zonen[1:176,1],col="red")
plot(pred22.julial, col="blue")
points(data2zonef1[1:176,1],col="red")
plot(pred21.juliam, col="blue",main="Medium")
points(data2zonen[1:41,81],col="red")
plot(pred22.juliam, col="blue")
points(data2zonef1[1:41,81],col="red")
plot(pred21.juliah, col="blue",main="High")
points(data2zonen[1:75,13],col="red")
plot(pred22.juliah, col="blue")
points(data2zonef1[1:75,13],col="red")
#plot all data
quartz()
par(mfrow=c(3,2))
plot(c1.hat.2zonel, col="blue",main="Low Near",ylab="Toluene Concentrations N",xlab="Time")
points(pred21.julial, col="red")
points(data2zonen[1:176,1],col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements N","State-by-state update",
                                                            "PMMH"),pch=15,cex=0.75,box.lwd = 0)

plot(c2.hat.2zonel, col="blue",main="Low Far",ylab="Toluene Concentrations F",xlab="Time")
points(pred22.julial, col="red")
points(data2zonef1[1:176,1],col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements F","State-by-state update",
                                                            "PMMH"),pch=15,cex=0.75,box.lwd = 0)


plot(data2zonen[1:75,13],col="green",main="Medium Near",ylab="Toluene Concentrations N",xlab="Time")
points(c1.hat.2zoneh, col="blue")
points(pred21.juliah, col="red")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements N","State-by-state update",
                                                            "PMMH"),pch=15,cex=0.75,box.lwd = 0)

plot(c2.hat.2zoneh, col="blue",main="Medium Far",ylab="Toluene Concentrations F",xlab="Time")
points(pred22.juliah, col="red")
points(data2zonef1[1:75,13],col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements F","State-by-state update",
                                                            "PMMH"),pch=15,cex=0.75,box.lwd = 0)

plot(c1.hat.2zonem, col="blue",main="High Near",ylab="2-butanone Concentrations N",xlab="Time")
points(pred21.juliam, col="red")
points(data2zonen[1:41,81],col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements N","State-by-state update",
                                                            "PMMH"),pch=15,cex=0.75,box.lwd = 0)

plot(c2.hat.2zonem, col="blue",main="High Far",ylab="2-butanone Concentrations F",xlab="Time")
points(pred22.juliam, col="red")
points(data2zonef1[1:41,81],col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements F","State-by-state update",
                                                            "PMMH"),pch=15,cex=0.75,box.lwd = 0)

#Turbulrnt eddy diffusion model
Dt=1
nloc<-5
ntime<-100
h=1
Gprime=351.4
x<-seq(1,nloc,length.out=ntime)
y<-seq(1,nloc,length.out=ntime)
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
    wt<-rgamma(ntime,aomega,bomega)
    csp[j,i]<-Gprime/(2*pi*Dt*(matdist[j]))*(1-
                                               (integrate(erf,lower=0,upper=(matdist[j])/sqrt(4*Dt*h*i)))$value)+wt[i]
  } 
}
Gprime/(2*pi*Dt*matdist[3])
#quartz()
plot(csp[1,])
points(csp[2,])
points(csp[3,])
points(csp[4,])
points(csp[5,])


csp2<-matrix(0,nloc,ntime)
csp2[,1]<-csp[,1]
for(i in 2:ntime){
  for(j in 1:nloc){
    set.seed(0000)
    wt<-rgamma(ntime,aomega,bomega)
    csp2[,i]<-csp2[,(i-1)]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
      (exp(-(matdist[j])^2/(4*Dt*i)))+(wt[i]-wt[(i-1)])
  }
}
csp2
quartz()
plot(csp[1,], csp2[1,])
abline(c(0,1))

#quartz()
plot(csp2[1,], ylim=c(0,max(csp2[1,])))
points(csp2[2,])
points(csp2[3,])
points(csp2[4,])
points(csp2[5,])

phi<-1
sigma<-1
tausq<-0.25
set.seed(123)
eta<-rnorm(ntime,0,0.1)
vt<-NULL
vt[1]<-0
for(i in 2:ntime){
  vt[i]<-vt[i-1]+eta[i]
}
vt
#library(geoR)
Sigma.exp <- function(x,y,phi,sigma) {
  Sigma <- matrix(rep(0, length(x)*length(y)), nrow=length(x))
  Sigma<- as.matrix(sigma*exp(-(dist(data.frame(cbind(x,y)),upper=TRUE, diag=TRUE,method="euclidean")*phi)))
  return(Sigma)
}
sigma.sp<-Sigma.exp(x,y,1,1)
dim(sigma.sp)
set.seed(123)
et<-rnorm(5,rep(0,5),c((tausq*diag(5)+sigma.sp[1:5,1:5])[,1]))
logyt3<-matrix(0,nrow=5, ncol=ntime)
for(i in 1:ntime){
  for(j in 1:5){
    logyt3[j,i]<-log(csp[j,i]+et[j])+eta[i]
  }
}

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
    csp[j,i+1]<-csp[j,i]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[j])^2/(4*Dt*i)))+(omega[i+1]-omega[i])
    logyt[j,i+1]~dnorm(log(csp[j,i+1]+et[j])+eta[i+1], tausq)
    }
    }
    for(i in 1:(N+1)){
    eta[i]~dnorm(0,10)
    #vt[(i+1)]<-vt[i]+eta[i]
    }
    for(j in 1:nloc){
    for(m in 1:nloc){
    k[m,j]<-sigma*exp(-phi*dist[m,j])
    kinv[m,j]<-inverse(k[m,j])
    }
    }
    for(j in 1:nloc){
    et[j]~dnorm(muet,kinv[j,1])
    #logyt[j,1]~dnorm(log(csp[j,1]+et[j])+eta[1], tausq)
    csp[j,1]<-exp(logyt[j,1])+omega[1]
    }
    omega[1]~dgamma(aomega,bomega)
    #vt[1]<-0
    Gprime~dunif(281,482)
    Dt~dunif(0,3)
    sigma~dgamma(2,2)
    phi~dunif(0.5,3)
    aomega~dunif(1,3)
    bomega~dunif(1/2,3)
    }",file="modeljagseddy.jag")
jags.eddy<-jags.model("modeljagseddy.jag",data=list("logyt"=logyt3, N=99,nloc=5,"h"=h,"muet"=0,
                                                    "tausq"=1/0.25,"pi"=pi, "dist"=distance, "matdist"=matdist))
update(jags.eddy, 1000)
mcmcjags.eddy<-jags.samples(jags.eddy,
                            c('Dt','Gprime','sigma','phi',"aomega","bomega",'csp'),
                            1000)
mcmcsamplesjags.eddy<-coda.samples(jags.eddy,
                                   c('Dt','Gprime','sigma','phi',"aomega","bomega",'csp'),
                                   1000)

quartz()
plot(mcmcsamplesjags.eddy[[1]][,1])
summary(mcmcsamplesjags.eddy)
m1.mcmc.eddy<-(as.mcmc(mcmcsamplesjags.eddy))
m1.mat.eddy<-as.matrix(mcmcsamplesjags.eddy)
m1.dat.eddy<-as.data.frame(m1.mat.eddy)
Dt.post.eddy<-m1.dat.eddy$Dt
gprime.post.eddy<-m1.dat.eddy$Gprime
aomega.post.eddy<-m1.dat.eddy$aomega
bomega.post.eddy<-m1.dat.eddy$bomega
sigma.post.eddy<-m1.dat.eddy$sigma
phi.post.eddy<-m1.dat.eddy$phi
quantile(bomega.post.eddy,c(0.025,0.25,0.5,0.75,0.975))

csp.hat<-matrix(0, nrow=5, ncol=100)
for(j in 1:5){
  for(i in 1:100){
    csp.hat[j,i] <- mean(mcmcsamplesjags.eddy[[1]][,4+(5*i+1-(6-j))])
  }
}
csp.hat
plot(csp[5,],csp.hat[5,])
abline(0,1)
quartz()
par(mfrow=c(3,2))
plot(csp[1,],col="red")
points(csp.hat[1,],col="blue")
points(exp(logyt3[1,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct1","chat1","yt1"),bg="white", lty=1,cex=0.8)
plot(csp[2,],col="red")
points(csp.hat[2,],col="blue")
points(exp(logyt3[2,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct2","chat2","yt2"),bg="white", lty=1,cex=0.8)
plot(csp[3,],col="red")
points(csp.hat[3,],col="blue")
points(exp(logyt3[3,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct3","chat3","yt3"),bg="white", lty=1,cex=0.8)
plot(csp[4,],col="red")
points(csp.hat[4,],col="blue")
points(exp(logyt3[4,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct4","chat4","yt4"),bg="white", lty=1,cex=0.8)
plot(csp[5,],col="red")
points(csp.hat[5,],col="blue")
points(exp(logyt3[5,]),col="green")
legend("bottomright",col=c("red","blue","green"),text.font=4,legend=
         c("Ct5","chat5","yt5"),bg="white", lty=1,cex=0.8)
#pf results
logyteddy<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/logyteddy.csv")
cspeddy<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/cspeddy.csv")
thetaeddyjulia<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmceddy.csv")
apply(thetaeddyjulia, 2, function(x){quantile(x)})
xjuliaeddy1<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmceddy1.csv")
xjuliaeddy2<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmceddy2.csv")
xjuliaeddy3<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmceddy3.csv")
xjuliaeddy4<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmceddy4.csv")
xjuliaeddy5<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmceddy5.csv")

prededdy1.julia<-apply(xjuliaeddy1, 2, mean)
prededdy2.julia<-apply(xjuliaeddy2, 2, mean)
prededdy3.julia<-apply(xjuliaeddy3, 2, mean)
prededdy4.julia<-apply(xjuliaeddy4, 2, mean)
prededdy5.julia<-apply(xjuliaeddy5, 2, mean)
#all sim
quartz()
par(mfrow=c(3,2))
plot(csp[1,],col="blue", ylab="Concentrations",xlab="Time", main="State-by-state update L1")
points(csp.hat[1,],col="red")
points(exp(logyt3[1,]),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)

plot(unlist(cspeddy[1,]),col="blue", ylab="Concentrations",xlab="Time", main="PMMH L1")
points(prededdy1.julia, col="red")
points(unlist(exp(logyteddy[1,])),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)
plot(csp[1,],col="blue", ylab="Concentrations",xlab="Time", main="State-by-state update L2")
points(csp.hat[1,],col="red")
points(exp(logyt3[1,]),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)
plot(unlist(cspeddy[2,]),col="blue", ylab="Concentrations",xlab="Time", main="PMMH L2")
points(prededdy2.julia, col="red")
points(unlist(exp(logyteddy[2,])),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)
plot(csp[3,],col="blue", ylab="Concentrations",xlab="Time", main="State-by-state update L3")
points(csp.hat[3,],col="red")
points(exp(logyt3[3,]),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)
plot(unlist(cspeddy[3,]),col="blue", ylab="Concentrations",xlab="Time", main="PMMH L3")
points(prededdy3.julia, col="red")
points(unlist(exp(logyteddy[3,])),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)
plot(csp[4,],col="blue", ylab="Concentrations",xlab="Time", main="State-by-state update L4")
points(csp.hat[4,],col="red")
points(exp(logyt3[4,]),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)

plot(unlist(cspedddy[4,]),col="blue", ylab="Concentrations",xlab="Time", main="PMMH L4")
points(prededdy4.julia, col="red")
points(unlist(exp(logyteddy[4,])),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)
plot(csp[5,],col="blue", ylab="Concentrations",xlab="Time", main="State-by-state update L5")
points(csp.hat[5,],col="red")
points(exp(logyt3[5,]),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)

plot(unlist(cspeddy[5,]),col="blue", ylab="Concentrations",xlab="Time", main="PMMH L5")
points(prededdy5.julia, col="red")
points(unlist(exp(logyteddy[5,])),col="green")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","True",
                                                            "Estimated"),pch=15,cex=0.75,box.lwd = 0)

#data
setwd("/Users/n_a_abdallah/Desktop/spatial/Project2/acetone/")
acetone <- list.files(pattern = ".csv$")
acetonedata<-lapply(acetone, read.csv)
plot(acetonedata[[2]][,2])

#k=0.00484
matdist<-c(0.41,1.07)
distance<-matrix(c(0,(1.07-0.41),(1.07-0.41),0),nrow=2, ncol=2, byrow=TRUE)
h=1
acetone18<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/acetone/acetone18.csv")
N<-length(acetone18[,2])-1
logacetone00484<-matrix(0,N+1,2)
logacetone00484<-log(acetone18[,3:4]*(58)/(24.45)+1)
cat("model{
    for (i in 1:N){
    omega[i+1]~dgamma(aomega,bomega)
    for(j in 1:2){
    csp[i+1,j]<-csp[i,j]+h*(Gprime/(4*(pi*Dt*i)^1.5))*
    (exp(-(matdist[j])^2/(4*Dt*i)))+(omega[i+1]-omega[i])
    logacetone00484[i+1,j]~dnorm(log(csp[i+1,j]+et[j])+eta[i+1], tausq)
    }
    }
    for(i in 1:(N+1)){
    eta[i]~dnorm(0,10)
    #vt[(i+1)]<-vt[i]+eta[i]
    }
    for(j in 1:nloc){
    for(m in 1:nloc){
    k[m,j]<-sigma*exp(-phi*dist[m,j])
    kinv[m,j]<-inverse(k[m,j])
    }
    }
    for(j in 1:2){
    et[j]~dnorm(muet,kinv[j,1])
    #logyt[j,1]~dnorm(log(csp[1,j]+et[j])+eta[1], tausq)
    csp[1,j]<-exp(logacetone00484[1,j])+omega[1]
    }
    omega[1]~dgamma(aomega,bomega)
    #vt[1]<-0
    Gprime~dunif(100,200000)
    Dt~dunif(0,1)
    sigma~dgamma(2,2)
    phi~dunif(0.5,3)
    aomega~dunif(1,3)
    bomega~dunif(1/2,3)
    }",file="modeljagseddydata.jag")
jags.eddy<-jags.model("modeljagseddydata.jag",data=list("logacetone00484"=logacetone00484, N=N,nloc=2,"h"=h,"muet"=0,
                                                        "tausq"=1/0.25,"pi"=pi, "dist"=distance, "matdist"=matdist))
update(jags.eddy, 1000)
mcmcjags.eddy<-jags.samples(jags.eddy,
                            c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                            1000)
mcmcsamplesjags.eddy<-coda.samples(jags.eddy,
                                   c('Dt','Gprime','sigma','phi',"aomega",'bomega','csp'),
                                   1000)

#quartz()
#plot(mcmcsamplesjags.eddy[[1]][,1])
#summary(mcmcsamplesjags.eddy)
m1.mcmc.eddy<-(as.mcmc(mcmcsamplesjags.eddy))
m1.mat.eddy<-as.matrix(mcmcsamplesjags.eddy)
m1.dat.eddy<-as.data.frame(m1.mat.eddy)
aomega.post.eddy<-m1.dat.eddy$aomega
bomega.post.eddy<-m1.dat.eddy$bomega
Dt.post.eddy<-m1.dat.eddy$Dt
gprime.post.eddy<-m1.dat.eddy$Gprime
sigma.post.eddy<-m1.dat.eddy$sigma
phi.post.eddy<-m1.dat.eddy$phi
quantile(Dt.post.eddy,c(0.025,0.25,0.5,0.75,0.975))/120
quantile(gprime.post.eddy,c(0.025,0.25,0.5,0.75,0.975))

csp.hat<-matrix(0, nrow=(N+1), ncol=2)
for(j in 1:2){
  for(i in 1:(N+1)){
    csp.hat[i,1] <- mean(mcmcsamplesjags.eddy[[1]][,4+i])
    csp.hat[i,2] <- mean(mcmcsamplesjags.eddy[[1]][,4+N+i])
    
  }
}
csp.hat
quartz()
plot(acetone18[,3]*54/24,csp.hat[,1])
abline(0,1)
quartz()
par(mfrow=c(1,2))
plot(acetone18[,3]*54/24,col="green")
points(csp.hat[,1],col="blue")
legend("bottomright",col=c("green","blue"),text.font=4,legend=
         c("measured","predicted"),bg="white", lty=1,cex=0.8)
plot(acetone18[,4]*54/24,col="green")
points(csp.hat[,2],col="blue")
legend("bottomright",col=c("green","blue"),text.font=4,legend=
         c("measured","predicted"),bg="white", lty=1,cex=0.8)
#PF
logyteddy10048<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmceddy1data2.csv")
logyteddy20048<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmceddy2data2.csv")
thetaeddyjulia0048<-read.csv("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmceddydata2.csv")
apply(thetaeddyjulia0048, 2, function(x){quantile(x)})[,1]/120
apply(thetaeddyjulia0048, 2, function(x){quantile(x)})[,2]
prededdy1.julia0048<-apply(logyteddy10048, 2, mean)
prededdy2.julia0048<-apply(logyteddy20048, 2, mean)
quartz()
par(mfrow=c(1,2))
plot(prededdy1.julia0048)
points(exp(logacetone00484[,1]),col="red")
plot(prededdy2.julia0048)
points(exp(logacetone00484[,2]),col="red")

#both plots
quartz()
par(mfrow=c(1,2))
plot(acetone18[,3]*54/24,col="green",ylab="Acetone Concentrations L1",xlab="Time", main="Location 1")
points(prededdy1.julia0048, col="red")
points(csp.hat[,1],col="blue")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","State-by-state update",
                                                            "PMMH"),pch=15,cex=0.75,box.lwd = 0)
plot(acetone18[,4]*54/24,col="green",ylab="Acetone Concentrations L2",xlab="Time", main="Location 2")
points(prededdy2.julia0048,col="red")
points(csp.hat[,2],col="blue")
legend("bottomright",col=c("green","blue", "red"), legend=c("Measurements","State-by-state update",
                                                            "PMMH"),pch=15,cex=0.75,box.lwd = 0)

