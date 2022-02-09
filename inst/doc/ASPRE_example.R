## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=TRUE----------------------------------------------------------------

# Set random seed
seed=463825
set.seed(seed)

# Libraries
library(mvtnorm)
library(matrixStats)
library(mle.tools)
library(OptHoldoutSize)

# Save plot to file, or not
save_plot=FALSE

# Force redo: set to TRUE to regenerate all datasets from scratch
force_redo=FALSE

#### ASPRE-related settings

# Total individuals in trial; all data
n_aspre_total=58794

# Population untreated PRE prevalence
pi_PRE = 1426/58974

# Maximum score sensitivity amongst highest 10%: assumed to be that of ASPRE
sens_max = (138+194)/2707  # = 0.122645 , from abstract of Rolnik 2017 Ultrasound in O&G

# Intervene with aspirin on this proportion of individuals
pi_intervention=0.1

# Aspirin reduces PRE risk by approximately this much
alpha=0.37
SE_alpha=0.09

# Candidate values for n
nval=round(seq(500,30000,length=100))

# Parameter calculation for N
N=400000; SE_N=1500

# Parameter calculation for k1
NICE_sensitivity=0.2
pi_1=NICE_sensitivity*(239/8875)/pi_intervention
pi_0=(1-NICE_sensitivity)*(239/8875)/(1-pi_intervention)
SE_pi_1=sqrt(pi_1*(1-pi_1)/(8875*0.1))
SE_pi_0=sqrt(pi_0*(1-pi_0)/(8875*0.9))
k1=pi_0*(1-pi_intervention) + pi_1*pi_intervention*alpha

# Standard error for k1
pi_1_s=rnorm(1000,mean=pi_1,sd=SE_pi_1)
pi_0_s=rnorm(1000,mean=pi_0,sd=SE_pi_0)
alpha_s=rnorm(1000,mean=alpha,sd=SE_alpha)
SE_k1=sd(pi_0_s*(1-pi_intervention) + pi_1_s*pi_intervention*alpha_s)





## ----echo=TRUE----------------------------------------------------------------
# Parameters of true ASPRE dataset
data(params_aspre)

# Simulate random dataset
X=sim_random_aspre(n_aspre_total,params=params_aspre)
X1=add_aspre_interactions(X)

# Risk will be monotonic to ASPRE risk, but we will transform to match
#  population prevalence of PE and sensitivity of ASPRE score.
risk0=aspre(X1)

# Find a linear transformation ax+b of lrisk such that population prevalence
#  and expected sensitivity match. Suppose P(Y_i=1)=score_i
# Expected sensitivity = E_{Y|scores}(sens)
#                      = (1/(pi_intervention*n_aspre_total))*E{sum_{i:score(i)>thresh} [Y_i]}
#                      = (1/5879)*sum_{i:score(i)>thresh} [(score(i)])
lrisk0=logistic(risk0)
f_ab=function(ab) {
  a=ab[1]; b=ab[2]
  risk_ab=a*lrisk0 + b
  pop_prev=mean(logit(risk_ab))
  q_pi=quantile(risk_ab,0.9)
  sens=(1/(pi_intervention*n_aspre_total))*sum(logit(risk_ab)*(risk_ab>q_pi))
  return((pop_prev-pi_PRE)^2 + (sens - sens_max)^2)
}
abmin=optim(c(1,0),f_ab)$par
lrisk=abmin[1]*lrisk0 +abmin[2]
risk=logit(lrisk)

# PRE is a 0/1 variable indicating whether that simulated patient had PRE.
PRE=rbinom(n_aspre_total,1,prob=risk) # ASPRE=ground truth


## ----echo=T,eval=F------------------------------------------------------------
#  set.seed(487276)
#  
#  # Start with estimates of k2 at 10 values of n
#  nn_par=round(runif(20,20,150)^2)
#  k2_par=0*nn_par;
#  for (i in 1:length(nn_par)) {
#    k2_par[i]=aspre_k2(nn_par[i],X,PRE)
#  }
#  
#  # Candidate values for n
#  nval=round(seq(500,30000,length=100))
#  
#  # Starting value for theta
#  theta=powersolve_general(nn_par,k2_par)$par
#  theta_se=powersolve_se(nn_par,k2_par,init=theta)
#  
#  # Rough estimate for variance of k2
#  dvar0=var(k2_par-powerlaw(nn_par,theta))
#  s2_par=rep(dvar0,length(k2_par))
#  
#  ## Successively add new points
#  for (i in 1:100) {
#    nxn=next_n(nval,nn_par,k2_par,var_k2 = s2_par,N=N,k1=k1,nmed=10)
#    if (any(is.finite(nxn))) n_new=nval[which.min(nxn)] else n_new=sample(nval,1)
#    k2_new=aspre_k2(n_new,X,PRE)
#    nn_par=c(nn_par,n_new)
#    k2_par=c(k2_par,k2_new)
#    s2_par=c(s2_par,dvar0)
#    print(i)
#  }
#  
#  # Resample k2(n), to avoid double-dipping effect
#  for (i in 1:length(nn_par)) {
#    k2_par[i]=aspre_k2(nn_par[i],X,PRE)
#  }
#  
#  # Transform to total cost
#  cc_par=k1*nn_par + k2_par*(N-nn_par)
#  
#  # Save
#  aspre_parametric=list(nn_par=nn_par,k2_par=k2_par,s2_par=s2_par,cc_par=cc_par)
#  save(aspre_parametric,file="data/aspre_parametric.RData")
#  

## ----echo=T-------------------------------------------------------------------
# Load data
data(aspre_parametric)
for (i in 1:length(aspre_parametric)) assign(names(aspre_parametric)[i],aspre_parametric[[i]])

theta=powersolve_general(nn_par,k2_par)$par
theta_se=powersolve_se(nn_par,k2_par,init=theta)

print(theta)
print(theta_se)

## ----echo=T-------------------------------------------------------------------

# Optimal holdout set size and cost
optim_aspre=optimal_holdout_size(N,k1,theta)
OHS_ASPRE=optim_aspre$size
MIN_COST_ASPRE=optim_aspre$cost

# Errors
cov_par=matrix(0,5,5);
cov_par[1,1]=SE_N^2; cov_par[2,2]=SE_k1^2
cov_par[3:5,3:5]=theta_se
CI_OHS_ASPRE=ci_ohs(N,k1,theta,sigma=cov_par,mode = "asymptotic",grad_nstar=grad_nstar_powerlaw,alpha = 0.1)

print(round(OHS_ASPRE))
print(round(MIN_COST_ASPRE))
print(round(CI_OHS_ASPRE))


## ----echo=F,fig.width=6,fig.height=6------------------------------------------

plot(0,xlim=range(nn_par),ylim=range(cc_par),type="n",
     xlab="Training set size",
     ylab=expression(paste("Total. cost ", "(", "","",
                           phantom() %prop% phantom(), " sens.", ")", "")))
points(nn_par,cc_par,pch=16,cex=0.5)
lines(nval,k1*nval + powerlaw(nval,theta)*(N-nval))
e_min=min(CI_OHS_ASPRE); e_max=max(CI_OHS_ASPRE); c_min=min(cc_par); c_max=max(cc_par);
polygon(c(e_min,e_min,e_max,e_max),c(c_min,c_max,c_max,c_min),
        col=rgb(1,0,0,alpha=0.2),border=NA)
points(OHS_ASPRE,MIN_COST_ASPRE,pch=16,col="red")

legend("topright",
       c("Cost function",
         "Est cost (d)",
         "OHS",
         "OHS err."),
       lty=c(1,NA,NA,NA),lwd=c(1,NA,NA,NA),pch=c(NA,16,16,16),pt.cex=c(NA,0.5,1,2),
       col=c("black","black","red",rgb(1,0,0,alpha=0.2)),bg="white",bty="n")

if (save_plot) dev.off()


## ----echo=T,eval=F------------------------------------------------------------
#  
#  plot(0,xlim=range(nn_par),ylim=range(cc_par),type="n",
#       xlab="Training set size",
#       ylab=expression(paste("Total. cost ", "(", "","",
#                             phantom() %prop% phantom(), " sens.", ")", "")))
#  points(nn_par,cc_par,pch=16,cex=0.5)
#  lines(nval,k1*nval + powerlaw(nval,theta)*(N-nval))
#  e_min=min(CI_OHS_ASPRE); e_max=max(CI_OHS_ASPRE); c_min=min(cc_par); c_max=max(cc_par);
#  polygon(c(e_min,e_min,e_max,e_max),c(c_min,c_max,c_max,c_min),
#          col=rgb(1,0,0,alpha=0.2),border=NA)
#  points(OHS_ASPRE,MIN_COST_ASPRE,pch=16,col="red")
#  
#  legend("topright",
#         c("Cost function",
#           "Est cost (d)",
#           "OHS",
#           "OHS err."),
#         lty=c(1,NA,NA,NA),lwd=c(1,NA,NA,NA),pch=c(NA,16,16,16),pt.cex=c(NA,0.5,1,2),
#         col=c("black","black","red",rgb(1,0,0,alpha=0.2)),bg="white",border=NA)
#  
#  if (save_plot) dev.off()
#  

## ----echo=T-------------------------------------------------------------------
# Variance and covariance parameters
var_u=1000
k_width=5000

## ----echo=T,eval=F------------------------------------------------------------
#  
#    # Begin as for parametric approach
#    set.seed(487276)
#  
#    # Start with estimates of k2 at 10 values of n
#    nn_emul=round(runif(20,20,150)^2)
#    k2_emul=0*nn_emul;
#    for (i in 1:length(nn_emul)) {
#      k2_emul[i]=aspre_k2(nn_emul[i],X,PRE)
#    }
#  
#    # Candidate values for n
#    nval=round(seq(500,30000,length=100))
#  
#    # Starting value for theta
#    theta=powersolve_general(nn_emul,k2_emul)$par
#  
#    # Rough estimate for variance of k2
#    dvar0=var(k2_emul-powerlaw(nn_emul,theta))
#    s2_emul=rep(dvar0,length(k2_emul))
#  
#  
#    ## Successively add new points
#    for (i in 1:100) {
#      nxn = exp_imp_fn(nval,nset=nn_emul,k2=k2_emul,var_k2=s2_emul,
#                       N=N,k1=k1,theta=theta,var_u=var_u,k_width=k_width)
#      n_new = nval[which.max(nxn)]
#      k2_new=aspre_k2(n_new,X,PRE)
#      nn_emul=c(nn_emul,n_new)
#      k2_emul=c(k2_emul,k2_new)
#      s2_emul=c(s2_emul,dvar0)
#      theta=powersolve_general(nn_emul,k2_emul)$par
#      print(c(i,n_new))
#    }
#  
#    # Transform estimated k2 to costs
#    cc_emul=k1*nn_emul + k2_emul*(N-nn_emul)
#  
#    # Save
#    aspre_emulation=list(nn_emul=nn_emul,k2_emul=k2_emul,s2_emul=s2_emul,cc_emul=cc_emul)
#    save(aspre_emulation,file="data/aspre_emulation.RData")
#  

## ----echo=T-------------------------------------------------------------------
# Load data
data(aspre_emulation)
for (i in 1:length(aspre_emulation)) assign(names(aspre_emulation)[i],aspre_emulation[[i]])

# Mean and variance of emulator for cost function, parametric assumptions satisfied
p_mu=mu_fn(nval,nset=nn_emul,k2=k2_emul,var_k2 = s2_emul,
           theta=theta,N=N,k1=k1,var_u=var_u,k_width=k_width)
p_var=psi_fn(nval,nset=nn_emul,var_k2=s2_emul,N=N,var_u=var_u,k_width=k_width)

OHS_ASPRE=nval[which.min(p_mu)]
MIN_COST_ASPRE=min(p_mu)
OHS_ERR=error_ohs_emulation(nn_emul,k2_emul,var_k2=s2_emul,N=N,k1=k1,alpha=0.1,
                            var_u=var_u,k_width=k_width,theta=theta)

print(round(OHS_ASPRE))
print(round(MIN_COST_ASPRE))
print(round(range(OHS_ERR)))

## ----echo=F,fig.width=6,fig.height=6------------------------------------------
plot(0,xlim=range(nn_emul),ylim=range(cc_emul),type="n",
     xlab="Training set size",
     ylab=expression(paste("Total. cost ", "(", "","",
                           phantom() %prop% phantom(), " sens.", ")", "")))
points(nn_emul,cc_emul,pch=16,cex=0.5)
lines(nval,p_mu)
lines(nval,p_mu+3*sqrt(pmax(0,p_var)),col="blue")
lines(nval,p_mu-3*sqrt(pmax(0,p_var)),col="blue")
e_min=min(OHS_ERR); e_max=max(OHS_ERR); c_min=min(cc_emul); c_max=max(cc_emul);
polygon(c(e_min,e_min,e_max,e_max),c(c_min,c_max,c_max,c_min),
        col=rgb(1,0,0,alpha=0.2),border=NA)
points(OHS_ASPRE,MIN_COST_ASPRE,pch=16,col="red")

legend("topright",
       c(expression(mu(n)),
         expression(mu(n) %+-% 3*sqrt(psi(n))),
         "Est cost (d)",
         "OHS",
         "OHS err."),
       lty=c(1,1,NA,NA,NA),lwd=c(1,1,NA,NA,NA),pch=c(NA,NA,16,16,16),pt.cex=c(NA,NA,0.5,1,1),
       col=c("black","blue","black","red",rgb(1,0,0,alpha=0.2)),bg="white",bty="n")

if (save_plot) dev.off()


## ----echo=T,eval=FALSE--------------------------------------------------------
#  plot(0,xlim=range(nn_emul),ylim=range(cc_emul),type="n",
#       xlab="Training set size",
#       ylab=expression(paste("Total. cost ", "(", "","",
#                             phantom() %prop% phantom(), " sens.", ")", "")))
#  points(nn_emul,cc_emul,pch=16,cex=0.5)
#  lines(nval,p_mu)
#  lines(nval,p_mu+3*sqrt(pmax(0,p_var)),col="blue")
#  lines(nval,p_mu-3*sqrt(pmax(0,p_var)),col="blue")
#  e_min=min(OHS_ERR); e_max=max(OHS_ERR); c_min=min(cc_emul); c_max=max(cc_emul);
#  polygon(c(e_min,e_min,e_max,e_max),c(c_min,c_max,c_max,c_min),
#          col=rgb(1,0,0,alpha=0.2),border=NA)
#  points(OHS_ASPRE,MIN_COST_ASPRE,pch=16,col="red")
#  
#  legend("topright",
#         c(expression(mu(n)),
#           expression(mu(n) %+-% 3*sqrt(psi(n))),
#           "Est cost (d)",
#           "OHS",
#           "OHS err."),
#         lty=c(1,1,NA,NA,NA),lwd=c(1,1,NA,NA,NA),pch=c(NA,NA,16,16,16),pt.cex=c(NA,NA,0.5,1,1),
#         col=c("black","blue","black","red",rgb(1,0,0,alpha=0.2)),bg="white",border=NA)
#  
#  if (save_plot) dev.off()
#  

