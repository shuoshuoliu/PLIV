library(MASS)

# Functions that contain the algorithm and how to simulate the data

##################################################
###### loglikelihood #####################
Lglik=function(x, y, z, c_len, t_len, theta){
  n=length(x)
  
  alpha=theta[1:(c_len+2)]                           #length c_len+1
  beta =theta[(c_len+3):(c_len+t_len+4)]             #length t_len+1
  c    =theta[(c_len+t_len+5):(2*c_len+t_len+4)]     #length c_len
  t    =theta[(2*c_len+t_len+5):(2*c_len+2*t_len+4)] #length t_len
  rho  =theta[length(theta)-1]                       #length 1
  sig2 =theta[length(theta)]                         #length 1
  
  xt=matrix(NA,n,t_len)
  for (j in 1:t_len){
    xt[,j]=ifelse(x>t[j],x-t[j],0)
  }
  xt=cbind(1,xt,x)  #beta*(1,(x,t)^+)
  
  zc=matrix(NA,n,c_len)
  for (i in 1:c_len){
    zc[,i]=ifelse(z>c[i],z-c[i],0)
  }
  zc=cbind(1,zc,z)
  
  betax=xt%*%beta   #[n,(j+1)]*(j+1)
  alphaz=zc%*%alpha #[n,(k+1)]*(k+1)
  u = y - betax
  v = x - alphaz
  
  temp=(u^2-2*rho*u*v+v^2)/sig2
  J = n*log(sig2)+n/2*log(1-rho^2)+1/2/(1-rho^2)*sum(temp) #nlminb optimizes the minimum
  return(J) 
}

################ Gradients ####################
Gradient = function(x, y, z, c_len, t_len, theta, info = FALSE){
  n=length(x)
  
  alpha=theta[1:(c_len+2)]                           #length c_len+1
  beta =theta[(c_len+3):(c_len+t_len+4)]             #length t_len+1
  c    =theta[(c_len+t_len+5):(2*c_len+t_len+4)]     #length c_len
  t    =theta[(2*c_len+t_len+5):(2*c_len+2*t_len+4)] #length t_len
  rho  =theta[length(theta)-1]                       #length 1
  sig2 =theta[length(theta)] 
  
  zc=matrix(NA,n,c_len)
  zcInd=zc  # the indicator
  for (i in 1:c_len){
    zc[,i]=ifelse(z>c[i],z-c[i],0)
    zcInd[,i]=ifelse(z>c[i],1,0)  
  }
  zc=cbind(1,zc,z) 
  
  xt=matrix(NA,n,t_len)
  xtInd=xt  # the indicator
  for (j in 1:t_len){
    xt[,j]=ifelse(x>t[j],x-t[j],0)
    xtInd[,j]=ifelse(x>t[j],1,0) #n*j
  }
  xt=cbind(1,xt,x)
  
  betax=xt%*%beta
  alphaz=zc%*%alpha
  u = y - betax
  v = x - alphaz
  
  alphazInd=t(t(zcInd)*alpha[2:(c_len+1)])
  betaxInd=t(t(xtInd)*beta[2:(t_len+1)])
  temp=(u^2-2*rho*u*v+v^2)

  gg_alpha=zc*as.vector(v-rho*u)/sig2/(1-rho^2)
  gg_beta=xt*as.vector(u-rho*v)/sig2/(1-rho^2)
  gg_c=alphazInd*as.vector(rho*u-v)/sig2/(1-rho^2)
  gg_t=betaxInd*as.vector(rho*v-u)/sig2/(1-rho^2)
  gg_rho=rho/(1-rho^2)-rho*temp/sig2/(1-rho^2)^2+u*v/sig2/(1-rho^2)
  gg_sig2=temp/2/sig2^2/(1-rho^2)-1/sig2
  grd=-c(colSums(gg_alpha),colSums(gg_beta),colSums(gg_c),colSums(gg_t),sum(gg_rho),sum(gg_sig2))
  
  if(info) {
    U=cbind(gg_alpha,gg_beta,gg_c,gg_t,gg_rho,gg_sig2)
    mtm = function(m) {return(m%*%t(m))}
    UU = apply(U, 1, mtm) ### p*p*n each with u_i*t(u_i)
    p=ncol(U)
    U2 = aperm(array(UU, c(p, p, n)), c(3, 1, 2))  ### now n*p*p
    pi_theta = colSums(U2, dims = 1)
    sd=sqrt(abs(diag(ginv(pi_theta))))

    ## used for robust only
    # for the alpha line
    #alpha_com=cbind(-zc,rho*xt,-rho*betaxInd,(2*rho*v-u-u*rho^2)/(1-rho^2),(rho*u-v)/sig2)
    #alpha_line=t(zc)%*%zc/sig2/(1-rho^2)
    #V_beta1 = mean(((u/sigma_u^2-rho*v/(sigma_u*sigma_v))*x1/(1-rho^2))^2)
    #print(sqrt(1/V_beta1))
    #var_beta1 = V_beta1
    
    #sd_beta1 = sqrt(1/V_beta1)
    #sd_beta1 = sqrt(V_beta1)
    #print(sd_beta1)
    ######### var c1
    #aa = rho*v/(sigma_u*sigma_v)-u/sigma_u^2
    
    #den1=density(x)
    #pt1=which(den1$x>=c1)[1]
    #f.x.c1 = 0.5*(den1$y[pt1]+den1$y[pt1-1])
    
    #cdf.x=ecdf(x)
    #F.x.c1=cdf.x(c1)
    
    #V_c1 = mean((beta1*aa*x11/(1-rho^2))^2)
    #var_c1 = -1/V_c1
    #pc1 = sum(beta1*(rho*v/(sigma_u*sigma_v)-u/sigma_u^2)*(1-F.x.c1)/(1-rho^2))
    #print(pc1)
    #var_c1 = t(pc1)%*%pc1
    #sd_c1 = sqrt(1/V_c1)
    #print(sd_c1)
    ######### var c2
    #bb = rho*u/(sigma_u*sigma_v)-v/sigma_v^2
    
    #den2=density(z)
    #pt2=which(den2$x>=c2)[1]
    #f.z.c2 = 0.5*(den2$y[pt2]+den2$y[pt2-1])
    
    #cdf.z=ecdf(z)
    #F.z.c2=cdf.z(c2)
    
    #V_c2 = mean((alpha1*bb*z11/(1-rho^2))^2)
    #var_c2 = -1/V_c2
    #sd_c2 = sqrt(1/V_c2)
    #print(sd_c2)
  return(sd)
  } else {
    return(grd)
  }
}

####################################
###gradient based method
####################################
GradFit = function(x, y, z, c_len, t_len) {
  
  theta0=TSRI(x, y, z, c_len, t_len)
  n=length(x)
  rho=theta0[length(theta0)-1]
  sig2=theta0[length(theta0)]
  
  #initial=c(alpha,beta,theta0$c_final,theta0$t_final,rho,sig2)
  
  lenlen=(c_len+t_len+4)
  bound=2
  lmax <- nlminb(theta0, Lglik, Gradient,
                 x=x, y=y, z=z, c_len=c_len, t_len=t_len,
                 lower=c(rep(-bound,c_len+2),rep(-bound,t_len+2),rep(min(z),c_len),
                         rep(min(x),t_len),-0.9999,0.0001),
                 upper=c(rep(bound,c_len+2),rep(bound,t_len+2),rep(max(z),c_len),
                         rep(max(x),t_len),0.9999,5))
  
    theta=lmax$par
    sd=Gradient(x, y, z, c_len, t_len, theta, info = TRUE)
    zscore = qnorm(1-.05/2)
    qtl = cbind(theta-zscore*sd, theta+zscore*sd)
    TAB1 = cbind(theta, sd, theta/sd, qtl[,1], qtl[,2], 2*pnorm(-abs(theta/sd)))
    colnames(TAB1) = c("Estimate", "Std.Err", "Z value", "95% CI(lower)", "95% CI(upper)", "Pr(>z)")
    #rownames(TAB1) = c(paste("alpha",0:(c_len+1)), paste("beta",0:(c_len+1)), 
                       #paste("c",1:c_len), paste("t",1:t_len), "rho", "sig2")
  
  result=list(theta = theta, conv=lmax$convergence,sd=sd,tab=TAB1)
  return(result)
}

#############################################################
###### two stage residual inclusion (this gives rough approx)
#############################################################
TSLS = function(x, y, z, c_len, t_len){
  n=length(x)
  cc = seq(min(z)+0.2, max(z)-0.2, by=0.1)
  t = seq(min(x)+0.2, max(x)-0.2, by=0.1)
  
  lenc=length(cc)
  lent=length(t)
  
  res_c=rep(NA,lenc)
  for (i in 1:lenc){
    zc=ifelse(z>cc[i],z-cc[i],0)
    res_c[i]=sum(residuals(lm(x ~ zc+z))^2)
  }
  c_final=cc[which.min(res_c)]
  
  res_t=rep(NA,lent)
  for (j in 1:lent){
    xt=ifelse(x>t[j],x-t[j],0)
    res_t[j]=sum(residuals(lm(y ~ xt+x))^2)
  }
  t_final=t[which.min(res_t)]
  
  ########
  zc=ifelse(z>c_final,z-c_final,0)
  xt=ifelse(x>t_final,x-t_final,0)

  s1=lm(x ~ zc+z)
  s2=lm(y ~ xt+x)
  alpha = as.vector(unlist(s1$coefficients))
  beta  = as.vector(unlist(s2$coefficients))#[-length(tsri2$coefficients)]
  res = residuals(s1)
  sig2=var(res)
  rho=cor(as.vector(x),res)
  
  alpha_se=as.vector(unlist(sqrt(diag(vcov(s1)))))
  beta_se=as.vector(unlist(sqrt(diag(vcov(s2)))))
  
  ans=list(alpha=alpha,beta=beta,c_final=c_final,t_final=t_final,
           rho=rho,sig2=sig2,alpha_se=alpha_se,beta_se=beta_se)
  return(ans)
}


TSRI = function(x, y, z, c_len, t_len){
  n=length(x)
  c = seq(min(z), max(z), length.out=c_len+2)
  c=c[-1]
  c=c[-length(c)]
  t = seq(min(x), max(x), length.out=t_len+2)
  t=t[-1]
  t=t[-length(t)]
  
  lenc=length(c)
  lent=length(t)
  
  zc=matrix(NA,n,lenc)
  for (i in 1:lenc){
    zc[,i]=ifelse(z>c[i],z-c[i],0)
  }
  zc=cbind(zc,z)
  xt=matrix(NA,n,lent)
  for (j in 1:lent){
    xt[,j]=ifelse(x>t[j],x-t[j],0)
  }
  
  tsri1 = lm(x ~ zc)
  res = residuals(tsri1)
  
  #xtt=cbind(xt,x,res)
  xtt=cbind(xt,x)
  tsri2 = lm(y ~ xtt)

  alpha = tsri1$coefficients
  beta  = tsri2$coefficients#[-length(tsri2$coefficients)]
  sig2=var(res)
  rho=cor(x,res)
  
  par=c(alpha,beta,c,t,rho,sig2)
  return(par)
}

#data generation function
data.gen =  function(n, alpha, beta, c, t, rho, sig2){
  sigma = matrix(c(sig2,sig2*rho,sig2*rho,sig2), ncol=2)
  error = mvtnorm::rmvnorm(n=n, mean=c(0,0), sigma=sigma)
  
  c_len=length(c)
  t_len=length(t)
  zc=matrix(NA,n,c_len)
  xt=zc
  
  z = rnorm(n, 0, 1)
  for (i in 1:c_len){
    zc[,i]=ifelse(z>c[i],z-c[i],0)
  }
  zc=cbind(1,zc,z)
  x=zc%*%alpha+error[,1]
  
  xt=matrix(NA,n,t_len)
  for (j in 1:t_len){
    xt[,j]=ifelse(x>t[j],x-t[j],0)
  }
  xt=cbind(1,xt,x)
  y = xt%*%beta + error[,2]
  
  result=list(y=y, x=x, z=z, zc=zc, xt=xt)
  return(result)
}