################################################################################
## R script for functions in OptHoldoutSize package                           ##
################################################################################
##
## Sami Haidar-Wehbe, Sam Emerson, James Liley
## October 2021
##

################################################################################
## Functions for ASPRE-related analyses                                       ##
################################################################################

##' Sensitivity at theshold quantile 10%
##'
##'
##' @export
##' @name sens10
##' @description Computes sensitivity of a risk score at a threshold at which 10% of samples (or some proportion pi_int) are above the threshold.
##' @keywords aspre
##' @param Y True labels (1 or 0)
##' @param Ypred Predictions (univariate; real numbers)
##' @param pi_int Compute sensitivity when a proportion pi_int of samples exceed threshold, default 0.1
##' @return Sensitivity at this threshold
##' @examples
##'
##' # Simulate
##' set.seed(32142)
##'
##' N=1000
##' X=rnorm(N); Y=rbinom(N,1,prob=logit(X/2))
##'
##' pi_int=0.1
##' q10=quantile(X,1-pi_int) # 10% of X values are above this threshold
##'
##' print(length(which(Y==1 & X>q10))/length(which(X>q10)))
##' print(sens10(Y,X,pi_int))
##'
sens10=function(Y,Ypred,pi_int=0.1) {
  qa=quantile(Ypred,1-pi_int)
  return(mean(Y*(Ypred >= qa))/pi_int)
}




##' Simulate random dataset similar to ASPRE training data
##'
##'
##' @export
##' @name sim_random_aspre
##' @description Generate random population of individuals (e.g., newly pregnant women) with given population parameters
##'
##' Assumes independence of parameter variation. This is not a realistic assumption, but is satisfactory for our purposes.
##' @keywords aspre
##' @param n size of population
##' @param params list of parameters
##' @return Matrix of samples
##' @examples
##'
##' # Load ASPRE related data
##' data(params_aspre)
##'
##' X=sim_random_aspre(1000,params_aspre)
##'
##' print(c(median(X$age),params_aspre$age$median))
##'
##' print(rbind(table(X$parity)/1000,params_aspre$parity$freq))
##'
sim_random_aspre=function(n,params=NULL) {
  if (is.null(params)) stop("Specify parameter params: this should be the object params_aspre from data(params_aspre)")
  out=list()
  ind=1
  for (i in 1:length(params)) {
    if (all(names(params[[i]])==c("median","IQR"))) {
      vals=rnorm(n,mean=params[[i]]$median,sd=(params[[i]]$IQR[2]-params[[i]]$median)/qnorm(0.75))
      out[[i]]=vals
    }
    if (all(names(params[[i]])==c("name","freq"))) {
      vals=factor(sample(params[[i]]$name,n,replace=TRUE,prob=params[[i]]$freq)) # note - strings not factors
      out[[i]]=vals
    }
    if (all(names(params[[i]])=="freq")) {
      vals=sample(0:1,n,replace=TRUE,prob=c(1-params[[i]]$freq,params[[i]]$freq))
      out[[i]]=vals
    }
  }
  out=as.data.frame(out)
  colnames(out)=names(params)
  return(out)
}


##' Add interaction terms corresponding to ASPRE model
##'
##'
##' @export
##' @name add_aspre_interactions
##' @keywords aspre
##' @description Add various interaction terms to X. Interaction terms correspond to those in ASPRE.
##' @param X data frame
##' @return New data frame containing interaction terms.
##' @examples
##'
##' # Load ASPRE related data
##' data(params_aspre)
##'
##' X=sim_random_aspre(1000,params_aspre)
##' Xnew=add_aspre_interactions(X)
##'
##' print(colnames(X))
##' print(colnames(Xnew))
##'
add_aspre_interactions=function(X) {

  # Pre-processing
  X$age_x=pmax(X$age-35,0) # age-35 if age>35, else 0
  X$height_x=X$height-164
  X$weight_ch=(X$weight-69)*(X$htn_hist==0) # If history of hypertension, do not count weight
  X$fam_hx_ch=(X$fam_hx)*(X$htn_hist==0) # If history of hypertension, do not count family history
  X$dm_hist_ch=(X$dm_hist)*(X$htn_hist==0) # If history of hypertension, do not count history of diabetes

  X$ga_prev_PE2=(X$last_delivery_gestational_age-24)^2 *(X$parity=="Par_prev_PE") # 0 if no previous PE
  X$years_no_prev_PE_1=((X$years_since_last_pregnancy)^(-1))*(X$parity=="Par_no_prev_PE")# etc
  X$years_no_prev_PE_sqrt=(pmax(X$years_since_last_pregnancy,1)^(-0.5))*(X$parity=="Par_no_prev_PE")
  X$ga_no_prev_PE2=(X$last_delivery_gestational_age-24)^2 *(X$parity=="Par_no_prev_PE") # 0 if no previous PE

  return(X)
}


##' Computes ASPRE score
##'
##'
##' @export
##' @name aspre
##' @description Computes ASPRE model prediction on a matrix `X` of covariates
##'
##' Full ASPRE model from https://www.nejm.org/doi/suppl/10.1056/NEJMoa1704559/suppl_file/nejmoa1704559_appendix.pdf
##'
##' Model is to predict gestational age at PE; that is, a higher score indicates a lower PE risk, so coefficients are negated for model to predict PE risk.
##' @keywords aspre
##' @param X matrix, assumed to be output of sim_random_aspre with parameter params=params_aspre and transformed using add_aspre_interactions
##' @return vector of scores.
##' @examples
##'
##' # Load ASPRE related data
##' data(params_aspre)
##'
##' X=sim_random_aspre(1000,params_aspre)
##' Xnew=add_aspre_interactions(X)
##'
##' aspre_score=aspre(Xnew)
##'
##' plot(density(aspre_score))
##'
aspre=function(X) {

  X1pred=data.frame(
    age=X$age_x,
    height=X$height_x,
    afc=(X$race=="Afro_Caribbean"),
    sa=(X$race=="S_Asian"),
    ch=X$htn_hist,
    sle_apl=X$sle_apl_hist,
    invitro=(X$conception=="In_vitro"),
    parous_prev_PE=(X$parity=="Par_prev_PE"),
    ga_prev_PE2=X$ga_prev_PE2,
    parous_no_prev_PE=(X$parity=="Par_no_prev_PE"),
    years_no_prev_PE_1=X$years_no_prev_PE_1,
    years_no_prev_PE_sqrt=X$years_no_prev_PE_sqrt,
    ga_no_prev_PE2=X$ga_no_prev_PE2,
    weight_ch=X$weight_ch,
    fam_hx_ch=X$fam_hx_ch,
    dm_hist_ch=X$dm_hist_ch
  )

  Xxpred=data.frame(
    map=pmax(X$map,0.1),
    uterine_artery_PI=pmax(X$uterine_artery_PI,0.1),
    serum_PPA=pmax(X$serum_PPA,0.1),
    serum_PPG=pmax(X$serum_PPG,0.1)
  )

  coefs=c(-0.206886, 0.11711, -2.6786, -1.129, -7.2897,
    -3.0519, -1.6327, -8.1667, 0.0271988, -4.335, -4.15137651,
    9.21473572, 0.01549673, -0.0694096, -1.7154, -3.3899
  )

  intercept_g=54.3637

  mu_g=(as.matrix(X1pred) %*% coefs) + intercept_g
  sd_g=6.8833

  beta0_x=c(0.09564,0.54453,-0.62165,-0.93687)
  beta1_x=c(-0.001824,-0.013143,0.014692,0.021930)


  cor_x=rbind(
    c(1,        -0.05133, -0.00497, -0.02791),
    c(-0.05133,        1, -0.15992, -0.15084),
    c(-0.00497, -0.15992,        1,  0.32085),
    c(-0.02791, -0.15084,  0.32085,        1))

  sd_x=c(0.03724,  0.12894,  0.23539,  0.17723)

  var_x=cor_x*outer(sd_x,sd_x)


  # Compute density at g
  dg=function(g) dmvnorm(log10(as.matrix(Xxpred)),mean=pmax(beta0_x+ beta1_x*g,0),sigma=var_x)*
    dnorm(g,mu_g,sd_g)

  # Integrate to find result
  numerator=rep(0,dim(X)[1]); denominator=rep(0,dim(X)[1])
  for (i in 24:37) numerator=numerator + dg(i)
  for (i in 38:100) denominator=denominator + dg(i)

  score=pmin(pmax(numerator/denominator,1e-3),0.95)

  # quick imputation for missing values, unimportant
  mscore=mean(score[which(!is.na(score))])
  score[which(is.na(score))]=mscore

  return(score)

}






##' Cost estimating function in ASPRE simulation
##'
##' @export
##' @name aspre_k2
##' @description Estimate cost at a given holdout set size in ASPRE model
##' @keywords aspre
##' @param n Holdout set size at which to estimate k_2 (cost)
##' @param X Matrix of predictors
##' @param PRE Vector indicating PRE incidence
##' @param seed Random seed; set before starting or set to NULL
##' @param pi_PRE Population prevalence of PRE if not prophylactically treated. Defaults to empirical value 1426/58974
##' @param pi_intervention Proportion of the population on which an intervention will be made. Defaults to 0.1
##' @param alpha Reduction in PRE risk with intervention. Defaults to empirical value 0.37
##' @return Estimated cost
##' @examples
##'
##' # Simulate
##' set.seed(32142)
##'
##' N=1000; p=15
##' X=matrix(rnorm(N*p),N,p); PRE=rbinom(N,1,prob=logit(X%*% rnorm(p)))
##' aspre_k2(1000,X,PRE)
aspre_k2=function(n,X,PRE,seed=NULL,pi_PRE=1426/58974,pi_intervention=0.1,alpha=0.37) {
  if (!is.null(seed)) set.seed(seed)
  n_aspre_total=length(PRE)

  sub=sample(n_aspre_total,n,replace=T) # subsample training data with replacement
  csub=setdiff(1:n_aspre_total,unique(sub)) # complement; use for testing performance

  # Subsets of X,Y
  Xsub=X[sub,]; PREsub=PRE[sub]
  Xc=X[csub,]; PREc=PRE[csub]

  # Manage missing factor levels
  for (i in 1:dim(Xc)[2]) {
    if ("factor" %in% class(Xc[,i])) {
      u=unique(Xsub[,i])
      w=which((Xc[,i] %in% u))
      Xc=Xc[w,]; PREc=PREc[w]
    }
  }

  # Predictor and predicted values
  g1=suppressWarnings(glm(PREsub~.,data=data.frame(Xsub,PREsub),family=binomial(link="logit")))
  score_p=suppressWarnings(predict(g1,data.frame(Xc),type="response"))

  # Evaluate cost
  dif=sens10(PREc,score_p)
  cost=pi_PRE - pi_intervention*alpha*dif

  return(cost)
}


