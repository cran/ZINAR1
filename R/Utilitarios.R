#' @importFrom gamlss.dist rPIG dPIG
#' @importFrom VGAM rzinegbin dzinegbin dzipois rzipois
#' @importFrom stats rbinom dbinom acf optim var
#' @importFrom scales percent
#' @importFrom grDevices dev.new
#' @import MASS
#' @import statmod
#' @import gtools
#' @import graphics



# ZI-INAR(1) process
########################

rzinar1=function(n,a,r,th,family)
{
  n_neg=100 #number of the initial values that will be de disregarded.
  muy=th[1]*(1-r)/(1-a)

  y=s=vector()
  y[1]=round(muy)

  if(family=="Po") w=rzipois(n+n_neg,th,r)
  if(family=="NB") w=rzinegbin(n+n_neg,pstr0=r,munb=th[1],size=th[2])
  if(family=="GI") w=rbinom(n+n_neg,1,1-r)*rPIG(n+n_neg,mu=th[1],sigma=1/th[2])

  for(i in 2:(n+n_neg))
  {
    s[i]=rbinom(1,y[i-1],a)
    y[i]=s[i]+w[i]
  }

  y=y[-(1:n_neg)]

  return(y)
}

# This function computes the moments estimation of the parameters
#################################################################

YW=function(y,model,family)
{

  alpha=max(0.01,acf(y,plot=F)$acf[2]) #autocorrelation of lag 1

  ################### estimation of rho

  delta=var(y)/mean(y)
  if(model=="zinar")
  {
    rho=0.5
    if(family=="Po")
    {
      rho=(delta-1)*(1+alpha)/( (delta-1)*(1+alpha) + mean(y)*(1-alpha) )
      rho=max(0.01,min(rho,0.99))
    }
  } else {
    rho=0
  }

  ################### estimation of the parameteres of the innovation.

  th=vector()
  th[1]=mean(y)*(1-alpha)/(1-rho) #estimate of mu

  if(family != "Po")
  {
    th[2]=th[1]/( (delta-1)*(1+alpha)-th[1]*rho ) #estimate of phi
    if(th[2]<=0 || th[2]>=5) th[2]=1
  }

  YW=c(alpha,rho,th)

  return(YW)

}

#ZIPIG DENSITY
##############

dpig=function(x,mu,sigma)
{
  phi=1/sigma

  #initially, we consider that all the values of the dpig function are zero.
  y=rep(0,length(x))

  #Then, we replace the zero values by the values of the ZIPIG density when x[i]>=0.
  if(any(x>=0))
  {
    pos=which(x>=0)
    y[pos]=dPIG(x[pos],mu=mu,sigma=1/phi)
  }

  return(y)
}

dzipig=function(x,mu,sigma,nu)
{
  phi=1/sigma

  #initially, we consider that all the values of the dzipig function are zero.
  y=rep(0,length(x))

  #Then, we replace the zero values by the value of the ZIPIG density when x[i]>=0.
  if(any(x>=0))
  {
    pos=which(x>=0)
    y[pos]=nu*(x[pos]==0)+(1-nu)*dPIG(x[pos],mu=mu,sigma=1/phi)
  }

  return(y)
}

# ZERO INFLATED DENSITY
#######################

dzi=function(z,r,th,family)
{

  if(family=="Po")
  {
    return( dzipois(z,th,r) )
  }

  if(family=="NB")
  {
    return( dzinegbin(z,pstr0=r,munb=th[1],size=th[2]) )
  }

  if(family=="GI")
  {
    return( dzipig(z,mu=th[1],sigma=1/th[2],nu=r) )
  }

}

#transition probability (ZI-INAR(1))
####################################

prob_trans=function(y,t,a,r,th,family)
{

  yn=y[t] #next state
  ylast=y[t-1] #current state

  top=min(yn,ylast)

  pr1=dbinom(0:top,ylast,a)
  pr2=dzi(yn-0:top,r,th,family)

  prob=sum(pr1*pr2)
  return(prob)
}

#LOG LIKELIHOOD FUNCTION (ZI-INAR(1))
#####################################

loglike=function(y,theta,family)
{
  n=length(y)

  alpha=theta[1]; rho=theta[2]; th=theta[3:length(theta)]

  L=sum(log(sapply(2:n,function(t){prob_trans(y,t,alpha,rho,th,family)}))) #log-likelihood

  return(L)
}


#distribution of (S[t]|y[t],y[t-1])
###################################

fst_propto=function(y,t,a,r,th,family)
{
  yn=y[t]
  ylast=y[(t-1)]

  top=min(yn,ylast)

  bin=dbinom(0:top,ylast,a)
  inov=dzi(yn-0:top,r,th,family)

  fx=bin*inov

  return(fx)
}

#distribution of (B[t]*S[t]|y[t],y[t-1])
########################################

fbtst_propto=function(y,t,a,r,th,family)
{
  yn=y[t]
  ylast=y[(t-1)]

  top=min(yn,ylast)

  bin=dbinom(0:top,ylast,a)
  inov=dzi(yn-0:top,0,th,family)

  fx=(1-r)*bin*inov

  return(fx)
}

#expectation
############

expectation=function(f)
{
  n=length(f)-1
  E=sum(0:n*f)
  return(E)
}

#LATENT (NEGATIVE BINOMIAL)
###########################

expectation_Gpropto=function(fbtst_p,y_obs,t,phi)
{
  fbtst_p=unlist(fbtst_p[[t]])
  top=length(fbtst_p)-1

  g=lgamma(y_obs[t]-0:top+phi)
  E=sum(g*fbtst_p)
  return(E)
}

#LOG LIKELIHOOD FUNCTION (ZINB-INAR(1))
#######################################

Qzinb=function(y_obs,fbtst_p,p_t,wt,btst,a,r,mu,phi)
{
  n=length(y_obs)+1

  ######### latent gt
  Eg_propto=function(t) expectation_Gpropto(fbtst_p,y_obs,t,phi)
  gt=sapply(1:(n-1), Eg_propto)
  gt=gt/p_t

  #Q function
  if(phi>0)
  {
    Q=sum( gt + (log(mu)-log(mu+phi))*((1-wt)*y_obs-btst)+
             ( phi*(log(phi)-log(mu+phi))-lgamma(phi) )*(1-wt) )
  } else {
    Q=-1e6
  }

  return(Q)
}

#LATENT (PIG)
#############

expectation_EZ=function(fst_p,y_obs,t,r,mu,phi)
{
  fst_p=unlist(fst_p[[t]])
  top=length(fst_p)-1
  st=0:top

  #expectation E(B[t]*Z[t]|s[t],y[t],y[t-1])
  Ez=(1-r)*(y_obs[t]-st+1)/mu*dpig(y_obs[t]-st+1,mu=mu,sigma=1/phi)/
    dzipig(y_obs[t]-st,mu=mu,sigma=1/phi,nu=r)

  vt=(y_obs[t]-st)+1*(y_obs[t]==st) #this is necessary to avoid 0/0 in the next expression

  #expectation E(B[t]/Z[t]|s[t],y[t],ty[t-1])
  Ezinv=(y_obs[t]>st)*mu/vt*dpig(y_obs[t]-st-1,mu,1/phi)/
    dpig(y_obs[t]-st,mu,1/phi)+
    (y_obs[t]==st)*(sqrt(phi*(phi+2*mu))+1)/phi*
    (1-r)*dpig(0,mu=mu,sigma=1/phi)/dzipig(0,mu=mu,sigma=1/phi,nu=r)

  E1=sum(Ez*fst_p) #numerator of the expectation E(B[t]*Z[t]|y[t],ty[t-1])
  E2=sum(Ezinv*fst_p) #numerator of the expectation E(B[t]/Z[t]|y[t],ty[t-1])

  return(c(E1,E2))
}

#ALGORITMO EM
#############

EM=function(y,initial,tol,int,model,family)
{

  n=length(y)

  L=vector() #vector of the log likelihood for the Aitken criteria.

  k=0 #count of the iterations

  a=initial[1]; r=ifelse(model=="inar",0,initial[2]); th=initial[3:length(initial)]
  theta_new=c(a,r,th)

  criterio=1
  while(criterio>tol && k<int)
  {
    k=k+1

    theta=theta_new
    if(k==1) L[1]=loglike(y,theta,family)
    L[k+1]=loglike(y,theta,family)

    p_t=sapply(2:n,function(t) prob_trans(y,t,a,r,th,family) ) #transition probability

    ############################ E STEP #############################
    #################################################################

    ############## latent st

    fst_p=lapply(2:n,function(t) fst_propto(y,t,a,r,th,family) )
    st=sapply(fst_p,expectation)
    st=st/p_t

    ############## latent btst

    fbtst_p=lapply(2:n,function(t) fbtst_propto(y,t,a,r,th,family) )
    btst=sapply(fbtst_p,expectation)
    btst=btst/p_t

    ############## latent wt

    if(model=="zinar")
    {
      wt=r*dbinom(y[2:n],y[1:(n-1)],a)
      wt=wt/p_t
    } else {
      wt=rep(0,n-1)
    }

    ############## latente bt*zt e bt/zt (gaussian inverse)

    if(family=="GI")
    {
      EZ_propto=function(t) expectation_EZ(fst_p,y[2:n],t,r,th[1],th[2])
      E=sapply(1:(n-1), EZ_propto)
      EZ=E[1,]/p_t
      EZinv=E[2,]/p_t
    }

    #################################################################
    ############################ M STEP #############################

    ######## innovation parameters

    if(family=="Po")
    {
      th=sum((1-wt)*y[2:n]-btst)/sum(1-wt)
    }

    if(family=="NB")
    {
      th[1]=mu=sum((1-wt)*y[2:n]-btst)/sum(1-wt)

      log_neg=function(phi) -Qzinb(y[2:n],fbtst_p,p_t,wt,btst,a,r,mu,phi)
      EMV=optim(par=th[2],fn=log_neg,method="BFGS")
      th[2]=EMV[[1]]

      th=c(th[[1]],th[[2]])
    }

    if(family=="GI")
    {
      th[1]=sum((1-wt)*y[2:n]-btst)/sum(EZ)
      th[2]=sum(1-wt)/(sum(EZ)+sum(EZinv)-2*sum(1-wt))
    }

    ######## alpha

    a=sum(st)/sum(y[1:(n-1)])

    ######## rho

    r=mean(wt)

    theta_new=c(a,r,th)

    ######################### Aitken criteria ##########################
    #############################################################

    L[k+2]=loglike(y,theta_new,family)
    ck=(L[k+2]-L[k+1])/(L[k+1]-L[k])
    L_inf=L[k+1]+(L[k+2]-L[k+1])/(1-ck)

    criterio=abs(L[k+2]-L_inf)

  }

  theta_new=rbind(round(theta_new,3))


  if(family=="Po")
  {
    if(model == "inar"){
      theta_new = rbind(theta_new[,c(1,3)])
      colnames(theta_new)=c("alpha","mu")
      rownames(theta_new)="EM estimates"
    } else {
      colnames(theta_new)=c("alpha","rho","mu")
      rownames(theta_new)="EM estimates"
    }

  } else {

    if(model == "inar"){
      theta_new = rbind(theta_new[,c(1,3,4)])
      colnames(theta_new)=c("alpha","mu","phi")
      rownames(theta_new)="EM estimates"
    } else {
      colnames(theta_new)=c("alpha","rho","mu","phi")
      rownames(theta_new)="EM estimates"
    }

  }

  return(list(parameters=theta_new,interactions=k))

}

#Exploratory
############

explore_zinar1 = function(x){

  if(sum(x%%1!=0)>0 | sum(x<0)>0){
    stop('x must be a ZINAR(1) process.')
  }

  df = data.frame(
    values = x
  )

  old.par = graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old.par))


  graphics::par(mfrow = c(2,2))

  stats::ts.plot(df$values,
                 main = 'Time Series Plot',
                 gpars = list(xlab = 'Time',
                              ylab = 'Value'))

  b = graphics::barplot(table(x)/length(x),
                        xlab = 'Values',
                        ylab = 'Percentage',
                        ylim = c(0,max(table(x)/length(x)) + 0.2),
                        main = 'Bar Chart',
                        cex.names = 0.6)
  text(b,table(x)/length(x),
       labels = percent(as.vector(table(x)/length(x))),
       pos = 3,
       srt = 90,
       offset = 1.5,
       cex = 0.9)


  stats::pacf(x, main = 'PACF')

  stats::acf(x, main = 'ACF')


}


