

```julia
Pkg.add("Rmath")
Pkg.build("Rmath")
Pkg.add("Distributions")
Pkg.add("RCall")
using Distributions
Pkg.add("Requests")
Pkg.add("DataFrames")
using Requests
using DataFrames
```


```julia
logyt=rand(Normal(),100)
n=100
Qprime=13.8
Gprime=351.54
V=100
sigma=0.5
aomega=2
bomega=1
kl=0.1
wt=rand(Gamma(aomega,bomega),n)

c=zeros(n)
for i in 1:(n-1)
  c[1]=1
  c[i+1]=(1-(Qprime+kl*V)/V)*c[i]+Gprime/V+wt[i+1]
end
srand(123)
vt=rand(Normal(0,sigma),n)
logyt=log(c)+(vt)
```




    100-element Array{Float64,1}:
     0.595134
     3.00345
     3.11458
     2.89186
     2.59154
     2.52275
     3.45231
     3.02925
     3.2489  
     2.96081
     2.91031
     2.64953
     2.66941
     ⋮       
     2.59575
     3.38919
     4.03901
     3.46275
     3.20035
     3.63186
     2.91012
     2.09172
     3.82615
     2.92951
     2.93079
     2.53731




```julia
function pf(nparticles::Integer,logyt::Array{Float64,1}, Qprime::Float64, Gprime::Float64, kl::Float64,aomega::Float64,bomega::Float64,sigma::Float64)
nparti=nparticles
cpartiold=zeros(nparti)
 cpartinew=zeros(nparti)
 wparti=zeros(nparti)
wpartin=zeros(nparti)
cpartires=zeros(nparti)
wpartires=zeros(nparti)
cestsir=zeros(100)
stdsir=zeros(100)
phat=zeros(100)
phat1=zeros(100)
mat=zeros(100)


nparti1=1/nparti
#distribution initial time
    cestsir[1]=logyt[1]
    for i in 1:nparti
          cpartiold[i]=exp(logyt[1])+rand(Gamma(1,2))
        end
      for j in 1:(n-1)
          for i in 1:nparti
    #prediction
            cpartinew[i]=(1-(Qprime+kl)/V)*cpartiold[i]+rand(Gamma(aomega,bomega))
    #weights using likelihood
            wparti[i]=exp(-((log(cpartinew[i])-logyt[j+1])/sigma)^2)
            end
      wtotal=sum(wparti)
      wpartin=wparti/wtotal
      phat[j]=(wtotal/nparti)
      #phat1=(phat1,phat)

  #update
      #cpartiw=cpartinew*wpartin
  #mean at time j+1
  #cestpf[j+1]<-cpartinew*wparti
  #resample
      cresa=zeros(nparti)
      uresa=zeros(nparti)
  #cumulative sum of weights
      cresa[1]=wpartin[1]
      for i in 2:nparti
        cresa[i]=cresa[i-1]+wpartin[i]
            end
      iresa=1
      uresa[1]=rand(Uniform(0,1))*nparti1
      for k in 1:nparti
        uresa[k]=uresa[1]+nparti1*(k-1)
        while uresa[k]>cresa[iresa]
          iresa=iresa+1
                end
        cpartires[k]=cpartinew[iresa]
        wpartires[k]=nparti1
                end
      cestsir[j+1]=mean(cpartires)
      stdsir[j+1]=sqrt(var(cpartires))
      cpartiold=cpartires
  #end of resampling
  #cpartiold<-cpartinew
    mat=cestsir
                    end
    return mat,phat
    end
V=100
@elapsed pzone1=pf(200,logyt,Qprime,Gprime,kl,aomega,bomega,sigma)
pzone1[1]
```

    WARNING: Method definition pf(Integer, Array{Float64, 1}, Float64, Float64, Float64, Float64, Float64, Float64) in module Main at In[45]:2 overwritten at In[215]:2.



    MethodError: no method matching pf(::Int64, ::Array{Float64,1}, ::Float64, ::Float64, ::Float64, ::Int64, ::Int64, ::Float64)
    Closest candidates are:
      pf(::Integer, ::Array{Float64,1}, ::Float64, ::Float64, ::Float64, ::Float64, ::Float64, ::Float64) at In[215]:2





```julia
niters=1000
tune=0.01
thin=10
theta=[Qprime, Gprime, kl, aomega,bomega,sigma]
p=length(theta)
thetamat=zeros(niters, p)
xmat=zeros(niters,100)
xtheta=zeros(100)
function pfmcmc1(niters, tune, thin,theta)
for i in 1:niters
  for j in 1:thin
    thetaprop=theta.*exp(rand(Normal(0,tune),p))
    priorprop=[pdf(Uniform(11,17),thetaprop[1]),pdf(Uniform(281,482),thetaprop[2]),pdf(Uniform(0,0.8),thetaprop[3]),
                     pdf(Uniform(0.5,3),thetaprop[4]),pdf(Uniform(0.5,3),thetaprop[5]),pdf(Gamma(1,2),thetaprop[6])]
    prior=[pdf(Uniform(11,17),theta[1]),pdf(Uniform(281,482),theta[2]),pdf(Uniform(0,0.8),theta[3]),
                     pdf(Uniform(0.5,3),theta[4]),pdf(Uniform(0.5,3),theta[5]),pdf(Gamma(1,2),theta[6])]
    pfprop=pf(200,logyt,thetaprop[1],thetaprop[2],thetaprop[3],thetaprop[4],thetaprop[5],thetaprop[6])
    llprop=log(prod(pfprop[2][1:99]))
    pftheta=pf(200,logyt,theta[1],theta[2],theta[3],theta[4],theta[5],theta[6])
    xtheta=pftheta[1]
    ll=log(prod(pftheta[2][1:99]))
    xprop=pfprop[1]
      if log(rand(Uniform()))< llprop+log(priorprop[1])+log(priorprop[2])+log(priorprop[3])+
         log(priorprop[4])+log(priorprop[5])+log(priorprop[6])-(ll+log(prior[1])+log(prior[2])+log(prior[3])+
                                                 log(prior[4])+log(prior[5])+log(prior[6]))
        theta=thetaprop
         ll=llprop
         xtheta=xprop
      else
      theta=theta
      ll=ll
      xtheta=xtheta
      end
  end
  thetamat[i,:]=theta
  xmat[i,:]=xtheta
end
return thetamat,xmat
end

```

    WARNING: Method definition pfmcmc1(Any, Any, Any, Any) in module Main at In[216]:9 overwritten at In[218]:10.





    pfmcmc1 (generic function with 1 method)




```julia
niters=100
@elapsed pfmcmcout=pfmcmc1(niters,tune, thin,theta)

```




23.13389615




```julia
thetamat=pfmcmcout[1]
xmat=pfmcmcout[2]
```


```julia
thetamatdat=DataFrame(thetamat)
xmatdat=DataFrame(xmat)

writetable("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmc.csv", thetamatdat)
writetable("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc.csv", xmatdat)


```


```julia
#two compartment model simulation
n=100
Qprime=13.8
Gprime=351.54
VF=1
VN=1
Sigma=diagm([1,1])
kl=0.1
beta=5
Var=diagm([1,1])
r=3
srand(123)
methods(MvNormal)
wt2=rand(MvLogNormal(2,1),n)
h=0.01
A=zeros(2,2)
A[1,1]=-(beta)/VN
A[1,2]= (beta)/VN
A[2,1]= beta/VF
A[2,2]= -(beta+Qprime)/VF+kl
g=zeros(1,2)
g[1,1]=Gprime/VN
c2zone=zeros(n,2)
for i in 1:(n-1)
  c2zone[1,:]=[0,0.5]
  c2zone[i+1,:]=(h*A+diagm([1,1]))*(c2zone[i,:])+(g*h)'+h*wt2[:,i+1]
end
c2zone

srand(123)
vt2=rand(MvNormal(2,1),n)
logyt2=log(c2zone)+(vt2)'
```


```julia
function pf2zone(nparticles,n,logyt, A,g,aomega,bomega,Sigma)
  #initialize vectors
  cestsir=zeros(n,2)
  cpartiold=zeros(nparticles,2)
  wparti=zeros(nparticles,1)
  cpartinew=zeros(nparticles,2)
  cpartires=zeros(nparticles,2)
  wpartires=zeros(nparticles,2)
  wpartin=zeros(nparticles,1)
  wtotal=wtotal2=cpartiw=cestpf=0.
  stdsir1=zeros(n,2)
  stdsir2=zeros(n,2)
  phat1=zeros(n,2)
  nparti=nparticles
  nparti1=1/nparti
  mat1=zeros(100)
  mat2=zeros(100)

  #distribution initial time
  cestsir[1,:]=exp(logyt2[1,:])
  for i in 1:nparti
    cpartiold[i,:]=exp(logyt2[1,:])+rand(Gamma(aomega,bomega),2)
  end
  #advance particles
  for j in 1:(n-1)
    for i in 1:nparti
      #prediction
      cpartinew[i,:]=(h*A+diagm([1,1]))*cpartiold[i,:]+(g*h)'+h*rand(Gamma(aomega,bomega),2)
      #weights using likelihood
      wparti[i,:]=exp(-((log(cpartinew[i,:])-logyt2[j+1,:])'*inv(Sigma)*(log(cpartinew[i,:])-logyt2[j+1,:]))/2)
    end
    wtotal=sum(wparti)
    wpartin=wparti/wtotal
    phat1[j]=wtotal/nparti
    #resample
    cresa=zeros(nparti)
    uresa=zeros(nparti)

    #cumulative sum of weights
    cresa[1]=wpartin[1]
    for i in 2:nparti
      cresa[i]=cresa[i-1]+wpartin[i]
    end
    iresa=1
    uresa[1]=rand(Uniform(0,1))*nparti1
    for k in 1:nparti
      uresa[k]=uresa[1]+nparti1*(k-1)
      while uresa[k]>cresa[iresa]
        iresa=iresa+1
        end
      cpartires[k,:]=cpartinew[iresa,:]
      wpartires[k,1]=nparti1

  end
    cestsir[j+1,1]=mean(cpartires[:,1])
    cestsir[j+1,2]=mean(cpartires[:,2])
    stdsir1[j+1]=sqrt(var(cpartires[:,1]))
    stdsir2[j+1]=sqrt(var(cpartires[:,2]))
    cpartiold=cpartires
    #end of resampling
    #cpartiold<-cpartinew
    mat1=(cestsir[:,1])
    mat2=(cestsir[:,2])
  end

  return mat1, mat2, phat1, phat2
  end

```

    WARNING: Method definition pf2zone(Any, Any, Any, Any, Any, Any, Any, Any) in module Main at In[93]:3 overwritten at In[95]:3.





    pf2zone (generic function with 1 method)




```julia
truepf2zone=pf2zone(200,100,logyt2, A,g,aomega,bomega,Sigma)
truepf2zone[1]
```




    0.356951296




```julia
function pfmcmc2(niters, tune, thin,theta)
for i in 1:niters
  for j in 1:thin
    thetaprop=theta.*exp(rand(Normal(0,tune),p))
    priorprop=[pdf(Uniform(11,17),thetaprop[1]),pdf(Uniform(281,482),thetaprop[2]),pdf(Uniform(0,0.8),thetaprop[3]),
                     pdf(Uniform(0.5,3),thetaprop[4]),pdf(Uniform(0.5,3),thetaprop[5]),pdf(Gamma(1,2),thetaprop[6][1,1]),
            pdf(Uniform(0,10),thetaprop[7])]
    prior=[pdf(Uniform(11,17),theta[1]),pdf(Uniform(281,482),theta[2]),pdf(Uniform(0,0.8),theta[3]),
                     pdf(Uniform(0.5,3),theta[4]),pdf(Uniform(0.5,3),theta[5]),pdf(Gamma(1,2),theta[6][1,1]),
            pdf(Uniform(0,10),theta[7])]
    A[1,1]=-(theta[7])/VN
    A[1,2]= (theta[7])/VN
    A[2,1]= theta[7]/VF
    A[2,2]= -(theta[7]+theta[1])/VF+theta[3]
    g[1,1]=theta[2]/VN
    pftheta=pf2zone(200,100,logyt2,A,g,theta[4],theta[5],theta[6])
    xtheta1=pftheta[1]
    xtheta2=pftheta[2]
    ll=log(prod(pftheta[3][1:99]))
    Aprop[1,1]=-(thetaprop[7])/VN
    Aprop[1,2]= (thetaprop[7])/VN
    Aprop[2,1]= thetaprop[7]/VF
    Aprop[2,2]= -(thetaprop[7]+thetaprop[1])/VF+thetaprop[3]
    gprop[1,1]=thetaprop[2]/VN
    pfthetaprop=pf2zone(200,100,logyt2,Aprop,gprop,thetaprop[4],thetaprop[5],thetaprop[6])
    xprop1=pfthetaprop[1]
    xprop2=pfthetaprop[2]
    llprop=log(prod(pfthetaprop[3][1:99]))
      if log(rand(Uniform()))< llprop+log(priorprop[1])+log(priorprop[2])+log(priorprop[3])+
         log(priorprop[4])+log(priorprop[5])+log(priorprop[6])+log(priorprop[7])-(ll+log(prior[1])+log(prior[2])+log(prior[3])+
                                                 log(prior[4])+log(prior[5])+log(prior[6])+log(prior[7]))
        theta=thetaprop
         ll=llprop
         xtheta1=xprop1
         xtheta2=xprop2
      else
          theta=theta
          ll=ll
          xtheta1=xtheta1
      xtheta2=xtheta2          
      end
  end
  thetamat[i,:]=theta
  xmat1[i,:]=xtheta1
  xmat2[i,:]=xtheta2
end
return thetamat,xmat1, xmat2
end

```

    WARNING: Method definition pfmcmc2(Any, Any, Any, Any) in module Main at In[190]:3 overwritten at In[220]:2.





    pfmcmc2 (generic function with 1 method)




```julia
 niters=1000
tune=0.01
thin=10
theta=[Qprime, Gprime, kl, aomega,bomega,Sigma[1,1], beta]
p=length(theta)
thetamat=zeros(niters, p)
xtheta1=zeros(niters, p)
xtheta2=zeros(niters, p)
xmat1=zeros(niters,100)
xmat2=zeros(niters,100)
A=zeros(2,2)
Aprop=zeros(2,2)
g=zeros(1,2)
gprop=zeros(1,2)


pfmcmcout2=pfmcmc2(niters,tune, thin,theta)


```




    (
    [13.3541 342.63 … 0.980302 5.13327; 13.9747 346.769 … 1.00587 4.79472; … ; 13.6578 433.395 … 19.689 6.56459; 13.8881 423.529 … 20.7916 6.31523],

    [0.0 4.09324 … 95.1632 95.3398; 0.0 4.09324 … 95.1632 95.3398; … ; 0.0 4.09324 … 95.1632 95.3398; 0.0 4.09324 … 95.1632 95.3398],

    [3.87689 4.07911 … 26.8296 26.8794; 3.87689 4.07911 … 26.8296 26.8794; … ; 3.87689 4.07911 … 26.8296 26.8794; 3.87689 4.07911 … 26.8296 26.8794])




```julia
for i in 1:niters
  for j in 1:thin
    thetaprop=theta.*exp(rand(Normal(0,tune),p))
    priorprop=[pdf(Uniform(11,17),thetaprop[1]),pdf(Uniform(281,482),thetaprop[2]),pdf(Uniform(0,0.8),thetaprop[3]),
                     pdf(Uniform(0.5,3),thetaprop[4]),pdf(Uniform(0.5,3),thetaprop[5]),pdf(Gamma(1,2),thetaprop[6][1,1]),
            pdf(Uniform(0,10),thetaprop[7])]
    prior=[pdf(Uniform(11,17),theta[1]),pdf(Uniform(281,482),theta[2]),pdf(Uniform(0,0.8),theta[3]),
                     pdf(Uniform(0.5,3),theta[4]),pdf(Uniform(0.5,3),theta[5]),pdf(Gamma(1,2),theta[6][1,1]),
            pdf(Uniform(0,10),theta[7])]
    A[1,1]=-(theta[7])/VN
    A[1,2]= (theta[7])/VN
    A[2,1]= theta[7]/VF
    A[2,2]= -(theta[7]+theta[1])/VF+theta[3]
    g[1,1]=theta[2]/VN
    pftheta=pf2zone(200,100,logyt2,A,g,theta[4],theta[5],theta[6])
    xtheta1=pftheta[1]
    xtheta2=pftheta[2]
    ll=log(prod(pftheta[3][1:99]))
    Aprop[1,1]=-(thetaprop[7])/VN
    Aprop[1,2]= (thetaprop[7])/VN
    Aprop[2,1]= thetaprop[7]/VF
    Aprop[2,2]= -(thetaprop[7]+thetaprop[1])/VF+thetaprop[3]
    gprop[1,1]=thetaprop[2]/VN
    pfthetaprop=pf2zone(200,100,logyt2,Aprop,gprop,thetaprop[4],thetaprop[5],thetaprop[6])
    xprop1=pfthetaprop[1]
    xprop2=pfthetaprop[2]
    llprop=log(prod(pfthetaprop[3][1:99]))
      if log(rand(Uniform()))< llprop+log(priorprop[1])+log(priorprop[2])+log(priorprop[3])+
         log(priorprop[4])+log(priorprop[5])+log(priorprop[6])+log(priorprop[7])-(ll+log(prior[1])+log(prior[2])+log(prior[3])+
                                                 log(prior[4])+log(prior[5])+log(prior[6])+log(prior[7]))
        theta=thetaprop
         ll=llprop
         xtheta1=xprop1
         xtheta2=xprop2
      else
          theta=theta
          ll=ll
          xtheta1=xtheta1
      xtheta2=xtheta2          
      end
  end
  thetamat[i,:]=theta
  xmat1[i,:]=xtheta1
  xmat2[i,:]=xtheta2
end
```


```julia
thetamatdat2=DataFrame(thetamat)
xmatdat1=DataFrame(xmat1)
xmatdat2=DataFrame(xmat2)
writetable("/Users/n_a_abdallah/Desktop/spatial/Project2/thetamcmc22.csv", thetamatdat2)
writetable("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc212.csv", xmatdat1)
writetable("/Users/n_a_abdallah/Desktop/spatial/Project2/xmcmc222.csv", xmatdat2)


```

    [1m[34mINFO: Nothing to be done
    [0m[1m[34mINFO: METADATA is out-of-date — you may not have the latest version of Requests
    [0m[1m[34mINFO: Use `Pkg.update()` to get the latest versions of your packages
    [0m[1m[34mINFO: Nothing to be done
    [0m[1m[34mINFO: METADATA is out-of-date — you may not have the latest version of DataFrames
    [0m[1m[34mINFO: Use `Pkg.update()` to get the latest versions of your packages
    [0mWARNING: Method definition describe(AbstractArray) in module StatsBase at /Users/n_a_abdallah/.julia/v0.5/StatsBase/src/scalarstats.jl:573 overwritten in module DataFrames at /Users/n_a_abdallah/.julia/v0.5/DataFrames/src/abstractdataframe/abstractdataframe.jl:407.



```julia

```
