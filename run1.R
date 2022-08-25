rm(list=ls())

source("main.R")
set.seed(200)
#generate true data#
n=500
alpha=c(-1,0.5,1)
beta=c(-0.2,0.5,0.5)
c=0.5
t=0
rho=0.5
s2=0.3
theta=c(alpha,beta,c,t,rho,s2)
data=data.gen(n, alpha, beta, c, t, rho, s2) #data generation
y=data$y
x=data$x
z=data$z
#summary(z)
#summary(x)

#####################################
B = 1000
out = matrix(0, B, 10)
out.sd = matrix(0, B, 10)

for (i in 1:B){ #repeat B times
  data=data.gen(n, alpha, beta, c, t, rho, s2)
  y=data$y
  x=data$x
  z=data$z
  
  fit=GradFit(x, y, z, 1, 1)
  out[i,] = fit$theta
  print(out[i,])
  out.sd[i,] = fit$sd
}

out.sd = out.sd[complete.cases(out.sd), ]
out = out[complete.cases(out), ]
nr=min(nrow(out.sd),nrow(out))
out=out[1:nr,]
out.sd=out.sd[1:nr,]

bias = colMeans(out) - theta
mc.sd = apply(out, 2, sd)
sd=colMeans(out.sd)
z=qnorm(0.975)

theta2=t(theta)[c(rep(1,nrow(out.sd))),]
cov=apply(((out-z*out.sd<=theta2) & (out+z*out.sd>=theta2)),2,sum)/nr
result=rbind(nr,bias,sd,mc.sd,cov)
result






