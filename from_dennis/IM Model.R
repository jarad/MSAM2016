library(plyr)
library(truncnorm)
library(ggplot2)
library(mnormt)
library(magic)
library(MCMCpack)
source("functions.R")

FG = read.csv("FGdata.csv",header = T)
g = read.csv("TakeIt14.csv",header = T)
FG = FG[FG$SEAS >= 2009,]
g = g[g$GID >= 2388,]
word = SmallN(FG,g,20)
FG = word$FG
g = word$g
g = g[g$pDiff >0 & g$DIST <= 76,]
FG$DIST = FG$DIST-18
g$DIST = g$DIST -18
g = arrange(g,GOOD)
o = glm(FG$GOOD~FG$FKICKER+FG$DIST*FG$FKICKER-1,family = "binomial")
oC = glm(FG$GOOD~FG$COND2+FG$SURF-1)
X = model.matrix(o)
Bp = dim(X)[2]/2
for(i in (dim(X)[2]/2+2):dim(X)[2]){X[,dim(X)[2]/2+1] = X[,dim(X)[2]/2+1]-X[,i]}
C = model.matrix(oC)
C = C[,-1]					
o2 = glm(g$TYPE~g$FKICKER+g$DIST*g$FKICKER-1)
oC2 = glm(g$TYPE~g$COND2+g$SURF-1)
X2 = model.matrix(o2)
for(i in (dim(X2)[2]/2+2):dim(X2)[2]){X2[,dim(X2)[2]/2+1] = X2[,dim(X2)[2]/2+1]-X2[,i]}
C2 = model.matrix(oC2)
C2 = C2[,-1]	
X2 = cbind(X2,C2)
X = cbind(X,C)
M = FG$GOOD
K = g$TYPE
D = g$pDiff



#Conjugate party
sims = 1100
pb <- txtProgressBar(min = 0, max = sims, style = 3)
#Ignore priors for now
#Starts and recording matrices
s1 = 1
s2 = 1
k1 = k2 = k3 = 0
r = 0
B = as.matrix(c(rnorm(Bp,14,1),rnorm(Bp,-0.14,0.1),0,0,0,0))
a = rep(1,4)
m = c(mean(B[1:(dim(X)[2]/2)]),mean(B[(dim(X)[2]/2+1):dim(X)[2]]))
Bdraws = matrix(rep(0,sims*dim(X)[2]),ncol = dim(X)[2])
mdraws = matrix(rep(0,sims*2),ncol = 2)
sdraws = matrix(rep(0,sims*2),ncol = 2)
rdraws = rep(0,sims)
adraws = matrix(rep(0,sims*4),ncol = 4)
Bdraws[1,] = B
mdraws[1,] = m
sdraws[1,1] = s1
sdraws[1,2] = s2
adraws[1,] = aK
Z = rep(0,length(K))
Y = rep(0,length(M))

i = 2
while(i<= sims)
{
  B = as.matrix(B,ncol = 1)
  #Hyper Prior updates
  #Mean term update
  #Sigt = Sigged(B[1:(2*Bp)],m)
  #sig = riwish(Bp+3,diag(2)+Sigt)
  #sigi = solve(sig)
  sig = matrix(c((1/s1)^2,r*(1/s1)*(1/s2),r*(1/s1)*(1/s2),(1/s2)^2),nrow = 2)
  sigi = solve(sig)
  V0 = diag(2)*100
  Vm = solve(solve(V0)+Bp*sigi)
  m = c(rmnorm(1,Vm%*%(Bp*sigi%*%c(mean(B[1:Bp]),mean(B[(Bp+1):(Bp*2)]))),Vm))
  
  #Variance terms updates
  #s1
  init = Hypelike(B,m,s1,s2,r)
  s1 = exp(rnorm(1,log(sdraws[i-1,1]),.5))
  finish = Hypelike(B,m,s1,s2,r)
  accept = exp(finish - init)
  if(accept < runif(1,0,1))
  {
    s1 = sdraws[i-1,1]
    k1 = k1+ 1
  }
  #s2
  init = Hypelike(B,m,s1,s2,r)
  s2 = exp(rnorm(1,log(sdraws[i-1,2]),.5))
  finish = Hypelike(B,m,s1,s2,r)
  accept = exp(finish - init)
  if(accept < runif(1,0,1))
  {
    s2 = sdraws[i-1,2]
    k2 = k2 + 1
  }
  
  #r
  init = Hypelike(B,m,s1,s2,r)
  r = rtruncnorm(1,a=-0.999,b=0.999,mean=rdraws[i-1],sd = 0.5)
  finish = Hypelike(B,m,s1,s2,r)
  accept = exp(finish - init)
  if(accept < runif(1,0,1))
  {
    r = rdraws[i-1]
    k3 = k3 + 1
  }
  
  #GIBBS B and Y
  ZK = (Z-a[1]-a[2]*D)/(a[3]+a[4]*D)
  XM = rbind(X,X2)
  XX = t(XM)%*%XM
  
  V = diag(c(rep((1/s1)^2,Bp),rep((1/s2)^2,Bp)))
  for(l in 1:(Bp)){
    V[l,Bp+l] = r*(1/s1)*(1/s2)
    V[Bp+l,l] = r*(1/s1)*(1/s2)
  }
  V = adiag(V,diag(4)*100)
  Vi = solve(V)
  B = as.matrix(B,ncol = 1)
  XB = X%*%B
  for(l in 1:length(M))
  {
    if(M[l] == 1){Y[l] = rtruncnorm(n=1,a = 0,b = Inf,mean=XB[l])}
    if(M[l] == 0){Y[l] = rtruncnorm(n=1,a = -Inf,b = 0,mean=XB[l])}
  }
  YZ = c(Y,ZK)
  #B = t(rmnorm(1,solve(XX)%*%(t(XM)%*%YZ),solve(XX)))
  B = c(rmnorm(1,solve(XX+Vi)%*%(Vi%*%c(rep(m[1],Bp),rep(m[2],Bp),0,0,0,0)+t(XM)%*%YZ),solve(XX+Vi)))
  B = as.matrix(B,ncol = 1)
  
  #GIBBS a
  Ya = X2%*%B
  Xa = as.matrix(data.frame("I" = rep(1,length(K)),"D" = D,"P" = Ya,"DP" = D*Ya))
  XXZ = solve(diag(4)/100+t(Xa)%*%Xa)
  XaB = Xa%*%a
  for(l in 1:length(K))
  {
    if(K[l] == 1){Z[l] = rtruncnorm(n=1,a = 0,b = Inf,mean=XaB[l])}
    if(K[l] == 0){Z[l] = rtruncnorm(n=1,a = -Inf,b = 0,mean=XaB[l])}
  }
  a = c(rmnorm(1,XXZ%*%(t(Xa)%*%Z),XXZ))
  Bdraws[i,] = B
  adraws[i,] = a
  sdraws[i,] = c(s1,s2)
  mdraws[i,] = m
  rdraws[i] = r
  setTxtProgressBar(pb, i)
  i = i + 1
}

#save.image("IM Draws.Rdata")

Beta = colMeans(Bdraws[100:(i-1),])
alpha = colMeans(adraws[100:(i-1),])




