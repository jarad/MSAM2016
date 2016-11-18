CreateX = function(data)
{
l = levels(data$FKICKER)
X = matrix(rep(0,length(l)),nrow = 1)
colnames(X) = l
p = data$player
for(i in 1:dim(data)[1])
{
	t = rep(0,length(l))
	j = 1
	while(p[i] != l[j]){j = j + 1}
	t[j] = 1
	X = rbind(X,t)
}
X = X[-1,]
X = X[,colSums(X) > 0]
return(X)
}

CreateXdist = function(data)
{
o = lm(data$GOOD~data$FKICKER)
l = levels(data$FKICKER)
X = matrix(rep(0,length(l)*2),nrow = 1)
colnames(X) = c(l,l)
p = data$FKICKER
for(i in 1:dim(data)[1])
{
	t = rep(0,length(l))
	j = 1
	while(p[i] != l[j]){j = j + 1}
	t[j] = 1
	t[j+length(l)] = FG$DIST[i]
	X = rbind(X,t)
}
X = X[-1,]
X = X[,colSums(X) > 0]
return(X)
}

SmallN = function(FG,g,n)
{
raw = ddply(FG,.(FKICKER),summarise, FG = mean(GOOD),number = length(GOOD))
enough = raw[raw$number < n,]
for(i in 1:dim(enough)[1])
{
FG = FG[FG$FKICKER != enough[i,1],]
g = g[g$FKICKER != enough[i,1],]
}
return(list("FG" = FG,"g" = g))
}

bern = function(Y,p)
{
temp = 0
for(i in 1:length(Y))
{
	this = 1/(1+exp(-p[i]))
	if(this == 1){out = 0.9999}
	if(this == 0){out = 0.0001}
	temp = temp + log((1-this)*(1-Y[i])+this*Y[i])
}
return(temp)
}

Binomial = function(x)
{
B = x[1:(dim(X)[2])]
g = x[(dim(X)[2]+1):(dim(X)[2]+dim(C)[2])]
a = x[(dim(X)[2]+dim(C)[2]+1):(dim(X)[2]+dim(C)[2]+4)]
psM = C%*%g+X%*%B
psK = C2%*%g+X2%*%B
out = bern(M,psM)+bern(K,a[1] + a[2]*D + a[3]*(1/(1+exp(-psK)))+a[4]*(D*(1/(1+exp(-psK)))))
return(-out)
}

stddevlike = function(B,m,rho,s)
{
muB = 0
sigB = 0
j = 1
while(j <= Bp)
{
	muB[j] = m[1]+s[1]/s[2]*rho*(B[j+Bp]-m[2])
	sigB[j] = (1-rho^2)*s[1]^2
	j = j + 1
}
while(j <= Bp*2)
{
	muB[j] = m[2]+s[2]/s[1]*rho*(B[j-Bp]-m[1])
	sigB[j] = (1-rho^2)*s[2]^2
	j = j + 1
}
out = dmnorm(B,muB,sigB*diag(Bp*2),log = TRUE)
return(out)
}

Sigged = function(B,m)
{
out = matrix(rep(0,4),nrow = 2)
for(j in 1:Bp)
{
	B = c(B)
	m = c(m)
	out= out +(B[c(j,j+Bp)]-m)%*%t(B[c(j,j+Bp)]-m)

}
return(out)
}


logit = function(x)
{
out = log(x/(1-x))
return(out)
}
rlogit = function(x)
{
out = 1/(1+exp(-x))
return(out)
}


sliceD = function(i,w1,w2,init)
{
finish = -Inf
u = log(runif(1,0,exp(init)))
U = runif(1,0,1)
L1 = B[i] - w1*U
R1 = L1 + w1
L2 = B[i+Bp] - w2*U
R2 = L2 + w2
while(fslice(i,L1) > u & fslice(i,R1) > u)
{
V = runif(1,0,1)
if(V < 0.5){L1 = L1 - (R1-L1)}
else{R1 = R1 + (R1-L1)}
}
while(fslice((i+Bp),L2) > u & fslice(i+Bp,R2) > u)
{
V = runif(1,0,1)
if(V < 0.5){L2 = L2 - (R2-L2)}
else{R2 = R2 + (R2-L2)}
}
while(finish < u)
{
B1 = runif(1,L1,R1)
B2 = runif(1,L2,R2)
Bt = B
Bt[i] = B1
Bt[i+Bp] = B2
psM = C%*%g+X%*%Bt
psK = C2%*%g+X2%*%Bt
finish = bern(M,psM)+bern(K,a[1] + a[2]*D + a[3]*(1/(1+exp(-psK)))+a[4]*(D*(1/(1+exp(-psK)))))+dmnorm(B[c(i,i+Bp)],m,sig,log = TRUE)
}
return(Bt)
}

fslice = function(i,B1)
{
Bt = B
Bt[i] = B1
psM = C%*%g+X%*%Bt
psK = C2%*%g+X2%*%Bt
if(i > Bp){d = dmnorm(B[c(i-Bp,i)],m,sig,log = TRUE)}
else{d = dmnorm(B[c(i,i+Bp)],m,sig,log = TRUE)}
finish = bern(M,psM)+bern(K,a[1] + a[2]*D + a[3]*(1/(1+exp(-psK)))+a[4]*(D*(1/(1+exp(-psK)))))+d
return(finish)
}

sliceD1 = function(i,w1,init)
{
finish = -Inf
u = log(runif(1,0,exp(init)))
U = runif(1,0,1)
L1 = B[i] - w1*U
R1 = L1 + w1
while(fslice(i,L1) > u & fslice(i,R1) > u)
{
V = runif(1,0,1)
if(V < 0.5){L1 = L1 - (R1-L1)}
else{R1 = R1 + (R1-L1)}
}
while(finish < u)
{
B1 = runif(1,L1,R1)
Bt = B
Bt[i] = B1
psM = C%*%g+X%*%Bt
psK = C2%*%g+X2%*%Bt
if(i > Bp){d = dmnorm(B[c(i-Bp,i)],m,sig,log = TRUE)}
else{d = dmnorm(B[c(i,i+Bp)],m,sig,log = TRUE)}
finish = bern(M,psM)+bern(K,a[1] + a[2]*D + a[3]*(1/(1+exp(-psK)))+a[4]*(D*(1/(1+exp(-psK)))))+d
}
return(Bt)
}

Hypelike = function(B,m,s1,s2,r)
{
sigi = pd.solve(matrix(c((1/s1)^2,r*(1/s1)*(1/s2),r*(1/s1)*(1/s2),(1/s2)^2),nrow = 2))
Vm = pd.solve(solve(V0)+Bp*sigi)
init = log(dmnorm(B[c(1,1+Bp)],m,Vm))
for(i in 2:Bp)
{
	init = init+log(dmnorm(B[c(i,i+Bp)],m,Vm))
}
return(init)
}

PG = function(z)
{
if(z<0.001){z = 0.001}
ok = 0
while(ok == 0)
{
out = 0
z = sqrt(z)
t = rnorm(1,0,1)^2
t = 1 + (t - sqrt(t*(4*z+t)))/(2*z)
ut = runif(1,0,1)
if(ut <= 1/(1+t))
{
	out = z/t
}
else{ out = z*t}
U = runif(1,0,1)
if(out > 4/3)
{
	ok = rightI(U,out)
}
else{	ok = leftI(U,out)}
}
return(out)
}


PG2 = function(z)
{
if(z<0.001){z = 0.001}
ok = 0
while(ok == 0)
{
out = rgig(n = 1, lambda = 0.5, chi = z,psi = 1)
U = runif(1,0,1)
if(out > 4/3)
{
	ok = rightI(U,out)
}
else{	ok = leftI(U,out)}
}
return(out)
}


rightI = function(U,out)
{
Z = 1
X = exp(-0.5*(out))
j = 0
while(j < 1000)
{
j = j + 1
Z = Z - (j+1)^2*X^((j+1)^2-1)
if(Z > U) {return(1)}
j = j + 1
Z = Z + (j+1)^2*X^((j+1)^2-1)
if(Z < U) {return(0)}
}
}

leftI = function(U,out)
{
H = log(2)/2 + 2.5*log(pi)-2.5*log(out)-(pi^2)/(2*out)+0.5*out
lu = log(U)
Z = 1
X = exp(-pi^2/(2*out))
K = out/pi^2
j = 0
while(j < 1000)
{
	j = j + 1
	Z = Z - K*X^(j^2-1)
	if(H+log(Z) > lu){return(1)}
	j = j + 1
	Z = Z + (j+1)^2*X^((j+1)^2-1)
	if(H+log(Z) < lu){return(0)}
}
}

