rm(list=ls())
library(ggplot2)
# needs sample size 2000 to have good coverage
# if sample small, then consider tight bound
source("main.R")
set.seed(200)
#generate true data#
n=500
alpha=c(-1,0.5,1)
beta=c(-0.2,1,0.5)
c=0.5
t=0
rho=0.2
s2=0.3
theta=c(alpha,beta,c,t,rho,s2)
data=data.gen(n, alpha, beta, c, t, rho, s2)
y=data$y
x=data$x
z=data$z

lx=ifelse(x>0,1,2)
lz=ifelse(z>0.5,1,2)

data=data.frame(x=x,y=y,z=z,lx=lx,lz=lz)

p1=ggplot(data, aes(x=z, y=x,color=lz)) + geom_point()+
  geom_smooth(method = "loess")+ geom_vline(xintercept=0.5, linetype="dashed", color = "red")+
  theme(legend.position = "none")
p2=ggplot(data, aes(x=x, y=y,color=lx)) + geom_point()+
  geom_smooth(method = "loess")+ geom_vline(xintercept=0, linetype="dashed", color = "red")+
  theme(legend.position = "none")

###########
set.seed(200)
#generate true data#
n=500
alpha=c(-1,0.5,1,1)
beta=c(-1,1.2,1,0.5)
c=c(-1,1)
t=c(-1,2)
rho=0.8
sig2=0.3

data=data.gen(n, alpha, beta, c, t, rho, sig2)
y=data$y
x=data$x
z=data$z

lx=ifelse(x>2,1,ifelse(x>-1&x<=2,2,3))
lz=ifelse(z>1,1,ifelse(z>-1&z<=1,2,3))

data=data.frame(x=x,y=y,z=z,lx=lx,lz=lz)

p3=ggplot(data, aes(x=z, y=x,color=lz)) + geom_point()+
  geom_smooth(method = "loess")+ geom_vline(xintercept=-1, linetype="dashed", color = "red")+
  geom_vline(xintercept=1, linetype="dashed", color = "red")+
  theme(legend.position = "none")
p4=ggplot(data, aes(x=x, y=y,color=lx)) + geom_point()+
  geom_smooth(method = "loess")+ geom_vline(xintercept=-1, linetype="dashed", color = "red")+
  geom_vline(xintercept=2, linetype="dashed", color = "red")+
  theme(legend.position = "none")

library(ggpubr)
library(grid)
figure=ggarrange(p1,p2,p3,p4, ncol = 2, nrow = 2)
#annotate_figure(figure, left = textGrob("MSPE", rot = 90, vjust = 1, gp = gpar(cex = 1.3)))

png("simu.png", width = 12, height = 8, units = 'in', res = 300)
figure
dev.off()

