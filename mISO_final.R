#######################################################################################################################
## Modified isotonic regression based phase I/II clinical trial designs identifying optimal biological dose
##          by Yingjie Qiu and Yong Zang
##
## This file contains six functions for implementation of the proposed designs to find the optimal 
## biological dose (OBD) for targted agents. The OBD is defined as the lowest dose with the
## highest rate of efficacy while safe.
##
## mISO.oc() and mISO_B.oc() are functions used to generate operating characteristics for the proposed designs
## mISO.df() and mISO_B.df() are functions used for dose-finding of the actual trial. 
## mISO.select() is a function used for determine the final OBD of the actual trial.
##
## At the end of this file, an example is provided to illustrate how to use the proposed designs to conduct a clinical trial.
########################################################################################################################

##########################################################################################################
## Function to generate operating characteristics of the modified isotonic design 
## To use this function, library "Iso" should be installed in R. 
##
##  Arguments:
## ttox : toxicity rates for each dose level
## teff: efficacy rate for each dose level
## cohortsize: sample size for each cohort
## ncohort: total number of cohort
## ntrial: the number of simulated trial
## phiT: upper bound of toxicity rate
## muT: threshold for posterior probability of toxicity; dose with toxicity probability larger than mut is excluded from the admissible set
## phiE: lower bound of efficacy rate
## muE: threshold for posterior probability of efficacy; dose with efficacy probability larger than muE is excluded from the admissible set
## alphaT: hyper parameters of prior of beta-binomial model of toxicity rate
## betaT: hyper parameters of prior of beta-binomial model of toxicity rate
## alphaE: hyper parameters of prior of beta-binomial model of efficacy rate
## betaE: hyper parameters of prior of beta-binomial model of efficacy rate
#########################################################################################################
mISO.oc = function(ttox,
                teff,
                cohortsize = 3,
                ncohort = 10,
                ntrial = 5000,
                phiT = 0.3,
                muT = 0.9,
                phiE = 0.5,
                muE = 0.85,
                alphaT=0.5,
                betaT=0.5,
                alphaE=0.5,
                betaE=0.5) {
  library("Iso")
  
  ###find admissible set of toxicity
  adm_tox <- function(n, ytox) {
    nn = n[n != 0]
    yytox = ytox[which(n != 0)]
    at = alphaT + yytox
    bt = betaT + nn - yytox
    Tox_prob = 1 - pbeta(phiT, at, bt)
    AT_naive = which(Tox_prob < muT)
    if (length(AT_naive)==0){
      AT=AT_naive
    }
    else{
      full_seq = seq(min(AT_naive),max(AT_naive),1)
      if (length(setdiff(full_seq, AT_naive)) == 0) {
        AT = AT_naive
      }
      else{
        AT = AT_naive[AT_naive < min(setdiff(full_seq, AT_naive))]
      }
    }
    return(AT)
  }
  
  ###find admissible set of efficacy
  adm_eff <- function(n, yeff) {
    nn = n[n != 0]
    yyeff = yeff[which(n != 0)]
    ae = alphaE + yyeff
    be = betaE + nn - yyeff
    Eff_prob = pbeta(phiE, ae, be)
    AR_naive = which(Eff_prob < muE)
    if (length(AR_naive)==0){
      AR=AR_naive
    }
    else{
      full_seq = seq(min(AR_naive),max(AR_naive),1)
      if (length(setdiff(full_seq, AR_naive)) == 0) {
        AR = AR_naive
      } else {
        AR = AR_naive[AR_naive > max(setdiff(full_seq, AR_naive))]
      }
    }
    return(AR)
  }
  
  
  ###AIC model selection
  
  AIC <- function(n, yeff) {
    aic = rep(0, length(n[n != 0]))
    for (l in 1:length(n[n != 0])) {
      if (l == 1) {
        jstar = length(n[n != 0])
        sm = sum(yeff[1:jstar])
        nm = sum(n[1:jstar])
        ql = sm / nm
        qfit = pava(y = ql, w = nm)
        ql_hat = qfit[1]
        likhood = (ql_hat ^ sm) * (1 - ql_hat) ^ (nm - sm)
        aic[l] = 2 * length(qfit) - 2 * log(likhood)
      } else{
        jstar = length(n[n != 0])
        sk = (yeff[which(n != 0)])[1:(l - 1)]
        nk = (n[n != 0])[1:(l - 1)]
        qk = sk / nk
        sm = sum(yeff[l:jstar])
        nm = sum(n[l:jstar])
        ql = sm / nm
        qfit = pava(y = c(qk, ql), w = c(nk, nm))
        qk_hat = qfit[1:(l - 1)]
        ql_hat = qfit[l]
        likhood = prod((qk_hat ^ sk) * (1 - qk_hat) ^ (nk - sk)) * (ql_hat ^
                                                                      sm) * (1 - ql_hat) ^ (nm - sm)
        aic[l] = 2 * length(qfit) - 2 * log(likhood)
      }
    }
    return(aic)
  }
  
  
  dselect = rep(0, ntrial)
  ndose = length(ttox)
  N = matrix(rep(0, ndose * ntrial), ncol = ndose)
  YTOX = matrix(rep(0, ndose * ntrial), ncol = ndose)
  YEFF = matrix(rep(0, ndose * ntrial), ncol = ndose)
  
  for (trial in 1:ntrial) {
    ytox = rep(0, ndose)
    yeff = rep(0, ndose)
    n = rep(0, ndose)
    d = 1 ##dose start from level 1
    
    ### dose-finding procedure
    for (i in 1:ncohort) {
      ytox[d] = ytox[d] + rbinom(1, cohortsize, ttox[d])
      yeff[d] = yeff[d] + rbinom(1, cohortsize, teff[d])
      n[d] = n[d] + cohortsize
      AT = adm_tox(n = n, ytox = ytox)
      AR = adm_eff(n = n, yeff = yeff)
      try = length(n[n != 0])
      if ((try %in% AT) & (try < ndose)) {
        d = d + 1
      } else {
        A = intersect(AT, AR)
        if (length(A) == 0) {
          d = 0
          break
        } else {
          OBD = A[which.min(AIC(n = n, yeff = yeff)[A])]
          if (OBD > d) {
            d = d + 1
          } else if (OBD < d) {
            d = d - 1
          } else {
            d = d
          }
        }
      }
    }
    if (d == 0) {
      dselect[trial] = 0
      N[trial, ] = n
      YTOX[trial, ] = ytox
      YEFF[trial, ] = yeff
    } else{
      AT = adm_tox(n = n, ytox = ytox)
      AR = adm_eff(n = n, yeff = yeff)
      A = intersect(AT, AR)
      if (length(A) == 0){
        dselect[trial] = 0
        N[trial, ] = n
        YTOX[trial, ] = ytox
        YEFF[trial, ] = yeff
      }
      else{
        dselect[trial] = A[which.min(AIC(n = n, yeff = yeff)[A])]
        N[trial, ] = n
        YTOX[trial, ] = ytox
        YEFF[trial, ] = yeff
      }
    }
  }
  selpercent = rep(0, ndose + 1)
  patpercent = matrix(rep(0, ntrial * ndose), ncol = ntrial, nrow = ndose)
  efficacy = rep(0, ntrial)
  toxicity = rep(0, ntrial)
  f <- function(x) {
    x[i] / sum(x)
  }
  ## Summarize results
  for (i in 0:ndose) {
    selpercent[(i + 1)] = sum(dselect == i) / ntrial * 100
  }
  print("selection probablity")
  cat(formatC(selpercent, digits = 1, format = "f"), sep = " ", "\n")
  for (i in 1:ndose) {
    patpercent[i, ] = apply(N, 1, f)
  }
  print("average percent of patients")
  cat(formatC(
    apply(patpercent, 1, mean) * 100,
    digits = 1,
    format = "f"
  ),
  sep = " ", "\n")
  print("average number of patients")
  cat(formatC(c(apply(N, 2, mean), sum(apply(
    N, 2, mean
  ))), digits = 1, format = "f"),
  sep = " ", "\n")
  print("average number of patients response to efficacy")
  cat(formatC(c(apply(YEFF, 2, mean), sum(
    apply(YEFF, 2, mean)
  )), digits = 1, format = "f"),
  sep = " ", "\n")
  for (i in 1:ntrial) {
    efficacy[i] = sum(YEFF[i, ]) / sum(N[i, ])
  }
  print("average percent of efficacy")
  cat(formatC(mean(efficacy) * 100, digits = 1, format = "f"),
      sep = " ", "\n")
  print("average number of patients response to toxicity")
  cat(formatC(c(apply(YTOX, 2, mean), sum(
    apply(YTOX, 2, mean)
  )), digits = 1, format = "f"),
  sep = " ", "\n")
  for (i in 1:ntrial) {
    toxicity[i] = sum(YTOX[i, ]) / sum(N[i, ])
  }
  print("average percent of toxicity")
  cat(formatC(mean(toxicity) * 100, digits = 1, format = "f"),
      sep = " ", "\n")
}

############################################################################
## Function for real trial dose-finding procedure of modified isotonic design
##
## n: number of patients treated at each dose level
## ytox: number of patients reporting toxicity at each dose level
## ytox: number of patients reporting efficacy at each dose level
## d: current dose level
############################################################################
mISO.df = function(n,ytox,yeff,d,phiT=0.3,muT=0.9,phiE=0.5,muE=0.85,
                   alphaT=0.5,betaT=0.5,alphaE=0.5,betaE=0.5){
  library("Iso")
  
  ###find admissible set of toxicity
  adm_tox <- function(n, ytox) {
    nn = n[n != 0]
    yytox = ytox[which(n != 0)]
    at = alphaT + yytox
    bt = betaT + nn - yytox
    Tox_prob = 1 - pbeta(phiT, at, bt)
    AT_naive = which(Tox_prob < muT)
    if (length(AT_naive)==0){
      AT=AT_naive
    }
    else{
      full_seq = seq(min(AT_naive),max(AT_naive),1)
      if (length(setdiff(full_seq, AT_naive)) == 0) {
        AT = AT_naive
      }
      else{
        AT = AT_naive[AT_naive < min(setdiff(full_seq, AT_naive))]
      }
    }
    return(AT)
  }
  
  ###find admissible set of efficacy
  adm_eff <- function(n, yeff) {
    nn = n[n != 0]
    yyeff = yeff[which(n != 0)]
    ae = alphaE + yyeff
    be = betaE + nn - yyeff
    Eff_prob = pbeta(phiE, ae, be)
    AR_naive = which(Eff_prob < muE)
    if (length(AR_naive)==0){
      AR=AR_naive
    }
    else{
      full_seq = seq(min(AR_naive),max(AR_naive),1)
      if (length(setdiff(full_seq, AR_naive)) == 0) {
        AR = AR_naive
      } else {
        AR = AR_naive[AR_naive > max(setdiff(full_seq, AR_naive))]
      }
    }
    return(AR)
  }
  
  
  ###AIC model selection
  
  AIC <- function(n, yeff) {
    aic = rep(0, length(n[n != 0]))
    for (l in 1:length(n[n != 0])) {
      if (l == 1) {
        jstar = length(n[n != 0])
        sm = sum(yeff[1:jstar])
        nm = sum(n[1:jstar])
        ql = sm / nm
        qfit = pava(y = ql, w = nm)
        ql_hat = qfit[1]
        likhood = (ql_hat ^ sm) * (1 - ql_hat) ^ (nm - sm)
        aic[l] = 2 * length(qfit) - 2 * log(likhood)
      } else{
        jstar = length(n[n != 0])
        sk = (yeff[which(n != 0)])[1:(l - 1)]
        nk = (n[n != 0])[1:(l - 1)]
        qk = sk / nk
        sm = sum(yeff[l:jstar])
        nm = sum(n[l:jstar])
        ql = sm / nm
        qfit = pava(y = c(qk, ql), w = c(nk, nm))
        qk_hat = qfit[1:(l - 1)]
        ql_hat = qfit[l]
        likhood = prod((qk_hat ^ sk) * (1 - qk_hat) ^ (nk - sk)) * (ql_hat ^
                                                                      sm) * (1 - ql_hat) ^ (nm - sm)
        aic[l] = 2 * length(qfit) - 2 * log(likhood)
      }
    }
    return(aic)
  }
  AT = adm_tox(n = n, ytox = ytox)
  AR = adm_eff(n = n, yeff = yeff)
  ndose = length(ytox)
  try = length(n[n != 0])
  if ((try %in% AT) & (try < ndose)) {
    d = d + 1
  } else {
    A = intersect(AT, AR)
    if (length(A) == 0) {
      d = 0
    } else {
      OBD = A[which.min(AIC(n = n, yeff = yeff)[A])]
      if (OBD > d) {
        d = d + 1
      } else if (OBD < d) {
        d = d - 1
      } else {
        d = d
      }
    }
  }
  return(list("next dose"=d))
}

############################################################################
## Function for real trial dose-recommendation of modified isotonic design
##
## If the trail does not stop early in the dose-finding procedure, 
##        we use this function to determine the final OBD
## n: number of patients treated at each dose level
## ytox: number of patients reporting toxicity at each dose level
## ytox: number of patients reporting efficacy at each dose level
############################################################################
mISO.select = function(n,ytox,yeff,phiT=0.3,muT=0.9,phiE=0.5,muE=0.85,
                   alphaT=0.5,betaT=0.5,alphaE=0.5,betaE=0.5){
  library("Iso")
  
  ###find admissible set of toxicity
  adm_tox <- function(n, ytox) {
    nn = n[n != 0]
    yytox = ytox[which(n != 0)]
    at = alphaT + yytox
    bt = betaT + nn - yytox
    Tox_prob = 1 - pbeta(phiT, at, bt)
    AT_naive = which(Tox_prob < muT)
    if (length(AT_naive)==0){
      AT=AT_naive
    }
    else{
      full_seq = seq(min(AT_naive),max(AT_naive),1)
      if (length(setdiff(full_seq, AT_naive)) == 0) {
        AT = AT_naive
      }
      else{
        AT = AT_naive[AT_naive < min(setdiff(full_seq, AT_naive))]
      }
    }
    return(AT)
  }
  
  ###find admissible set of efficacy
  adm_eff <- function(n, yeff) {
    nn = n[n != 0]
    yyeff = yeff[which(n != 0)]
    ae = alphaE + yyeff
    be = betaE + nn - yyeff
    Eff_prob = pbeta(phiE, ae, be)
    AR_naive = which(Eff_prob < muE)
    if (length(AR_naive)==0){
      AR=AR_naive
    }
    else{
      full_seq = seq(min(AR_naive),max(AR_naive),1)
      if (length(setdiff(full_seq, AR_naive)) == 0) {
        AR = AR_naive
      } else {
        AR = AR_naive[AR_naive > max(setdiff(full_seq, AR_naive))]
      }
    }
    return(AR)
  }
  
  
  ###AIC model selection
  
  AIC <- function(n, yeff) {
    aic = rep(0, length(n[n != 0]))
    for (l in 1:length(n[n != 0])) {
      if (l == 1) {
        jstar = length(n[n != 0])
        sm = sum(yeff[1:jstar])
        nm = sum(n[1:jstar])
        ql = sm / nm
        qfit = pava(y = ql, w = nm)
        ql_hat = qfit[1]
        likhood = (ql_hat ^ sm) * (1 - ql_hat) ^ (nm - sm)
        aic[l] = 2 * length(qfit) - 2 * log(likhood)
      } else{
        jstar = length(n[n != 0])
        sk = (yeff[which(n != 0)])[1:(l - 1)]
        nk = (n[n != 0])[1:(l - 1)]
        qk = sk / nk
        sm = sum(yeff[l:jstar])
        nm = sum(n[l:jstar])
        ql = sm / nm
        qfit = pava(y = c(qk, ql), w = c(nk, nm))
        qk_hat = qfit[1:(l - 1)]
        ql_hat = qfit[l]
        likhood = prod((qk_hat ^ sk) * (1 - qk_hat) ^ (nk - sk)) * (ql_hat ^
                                                                      sm) * (1 - ql_hat) ^ (nm - sm)
        aic[l] = 2 * length(qfit) - 2 * log(likhood)
      }
    }
    return(aic)
  }
  AT = adm_tox(n = n, ytox = ytox)
  AR = adm_eff(n = n, yeff = yeff)
  A = intersect(AT, AR)
  OBD = 0
  if (length(A) == 0){
    OBD = 0
  }
  else{
    OBD = A[which.min(AIC(n = n, yeff = yeff)[A])]
  }
  return(list("OBD"=OBD))
}

##########################################################################################################
## Function to generate operating characteristics of the mISO-B design 
## To use this function, library "Iso" should be installed in R. 
##
##  Arguments:
## ttox : toxicity rates for each dose level
## teff: efficacy rate for each dose level
## cohortsize: sample size for each cohort
## ncohort: total number of cohort
## ntrial: the number of simulated trial
## phiT: upper bound of toxicity rate
## muT: threshold for posterior probability of toxicity; dose with toxicity probability larger than mut is excluded from the admissible set
## phiE: lower bound of efficacy rate
## muE: threshold for posterior probability of efficacy; dose with efficacy probability larger than muE is excluded from the admissible set
## alphaT: hyper parameters of prior of beta-binomial model of toxicity rate
## betaT: hyper parameters of prior of beta-binomial model of toxicity rate
## alphaE: hyper parameters of prior of beta-binomial model of efficacy rate
## betaE: hyper parameters of prior of beta-binomial model of efficacy rate
## dist1:the underlying distribution of the time to efficacy outcomes; dist1=1 is the uniform distribution, 
##        dist1=2 corresponds to the Weibull distribution; dist1=3 is the log-logistic distribution.
## dist2: the underlying distribution of the time to toxicity outcomes; dist1=1 is the uniform distribution, 
##        dist1=2 corresponds to the Weibull distribution; dist1=3 is the log-logistic distribution.
## dist3: the underlying distribution of patient arrival time; dist3=1 is the uniform distribution, 
##        dist3=2 is the exponential distribution
## alpha1: a number from (0,1) that controls alpha*100% efficacy events in (0, 1/2T)
## alpha2: a number from (0,1) that controls alpha*100% toxicity events in (0, 1/2T)
## maxt1: time assess window for efficacy
## maxt2: time assess window for toxicity
## accrual: the accrual rate, i.e., the number of patients accrued in 1 unit of time
#########################################################################################################
mISO_B.oc <- function(ttox,
                   teff,
                   cohortsize = 3,
                   ncohort = 20,
                   ntrial = 5000,
                   phiT = 0.3,
                   muT = 0.95,
                   phiE = 0.5,
                   muE = 0.8,
                   alphaT=0.5,
                   betaT=0.5,
                   alphaE=0.5,
                   betaE=0.5,
                   dist1 = 2,# distribution of efficacy
                   dist2 = 2,#distribution of toxicity
                   dist3 = 1,#patient arrival
                   alpha1 = 0.5,
                   alpha2=0.5,
                   maxt1 = 3,# time assess window for efficacy
                   maxt2=3,# time assess window for toxicity
                   accrual = 2) {
  library("Iso")
  
  ###find admissible set of toxicity
  adm_tox <- function(n, ytox) {
    nn = n[n != 0]
    yytox = ytox[which(n != 0)]
    at = alphaT + yytox
    bt = betaT + nn - yytox
    Tox_prob = 1 - pbeta(phiT, at, bt)
    AT_naive = which(Tox_prob < muT)
    if (length(AT_naive) == 0) {
      AT = AT_naive
    }
    else{
      full_seq = seq(min(AT_naive), max(AT_naive), 1)
      if (length(setdiff(full_seq, AT_naive)) == 0) {
        AT = AT_naive
      }
      else{
        AT = AT_naive[AT_naive < min(setdiff(full_seq, AT_naive))]
      }
    }
    return(AT)
  }
  
  ###find admissible set of efficacy
  adm_eff <- function(n, yeff) {
    nn = n[n != 0]
    yyeff = yeff[which(n != 0)]
    ae = alphaE + yyeff
    be = betaE + nn - yyeff
    Eff_prob = pbeta(phiE, ae, be)
    AR_naive = which(Eff_prob < muE)
    if (length(AR_naive) == 0) {
      AR = AR_naive
    }
    else{
      full_seq = seq(min(AR_naive), max(AR_naive), 1)
      if (length(setdiff(full_seq, AR_naive)) == 0) {
        AR = AR_naive
      } else {
        AR = AR_naive[AR_naive > max(setdiff(full_seq, AR_naive))]
      }
    }
    return(AR)
  }
  ###AIC model selection
  
  AIC <- function(n, yeff) {
    aic = rep(0, length(n[n != 0]))
    for (l in 1:length(n[n != 0])) {
      if (l == 1) {
        jstar = length(n[n != 0])
        sm = sum(yeff[1:jstar])
        nm = sum(n[1:jstar])
        ql = sm / nm
        qfit = pava(y = ql, w = nm)
        ql_hat = qfit[1]
        likhood = (ql_hat ^ sm) * (1 - ql_hat) ^ (nm - sm)
        aic[l] = 2 * length(qfit) - 2 * log(likhood)
      } else{
        jstar = length(n[n != 0])
        sk = (yeff[which(n != 0)])[1:(l - 1)]
        nk = (n[n != 0])[1:(l - 1)]
        qk = sk / nk
        sm = sum(yeff[l:jstar])
        nm = sum(n[l:jstar])
        ql = sm / nm
        qfit = pava(y = c(qk, ql), w = c(nk, nm))
        qk_hat = qfit[1:(l - 1)]
        ql_hat = qfit[l]
        likhood = prod((qk_hat ^ sk) * (1 - qk_hat) ^ (nk - sk)) * (ql_hat ^
                                                                      sm) * (1 - ql_hat) ^ (nm - sm)
        aic[l] = 2 * length(qfit) - 2 * log(likhood)
      }
    }
    return(aic)
  }
  gen.tite <- function(dist = 1,
                       n,
                       pi,
                       alpha = 0.5,
                       Tobs = 1) {
    ############ subroutines ############
    weib <- function(n, pi, pihalft)
    {
      ## solve parameters for Weibull given pi=1-S(T) and phalft=1-S(T/2)
      alpha = log(log(1 - pi) / log(1 - pihalft)) / log(2)
      
      lambda = -log(1 - pi) / (Tobs ^ alpha)
      
      t = (-log(runif(n)) / lambda) ^ (1 / alpha)
      
      return(t)
      
    }
    
    llogit <- function(n, pi, pihalft)
    {
      ## solve parameters for log-logistic given pi=1-S(T) and phalft=1-S(T/2)
      alpha = log((1 / (1 - pi) - 1) / (1 / (1 - pihalft) - 1)) / log(2)
      
      lambda = (1 / (1 - pi) - 1) / (Tobs ^ alpha)
      
      t = ((1 / runif(n) - 1) / lambda) ^ (1 / alpha)
      
      return(t)
      
    }
    ############ end of subroutines ############
    
    
    tox = rep(0, n)
    
    t.tox = rep(0, n)
    
    
    #### uniform
    if (dist == 1) {
      # 50% event in (0, 1/2T)
      tox = rbinom(n, 1, pi)
      
      ntox.st = sum(tox)
      
      t.tox[tox == 0] = Tobs
      
      t.tox[tox == 1] = runif(ntox.st, 0, Tobs)
      
    }
    #### Weibull
    if (dist == 2)
    {
      pihalft = alpha * pi
      # alpha*100% event in (0, 1/2T)
      t.tox = weib(n, pi, pihalft)
      
      tox[t.tox <= Tobs] = 1
      
      ntox.st = sum(tox)
      
      t.tox[tox == 0] = Tobs
      
    }
    #### log-logistic
    if (dist == 3)
    {
      pihalft = alpha * pi
      # alpha*100% event in (0, 1/2T)
      t.tox = llogit(n, pi, pihalft)
      
      tox[t.tox <= Tobs] = 1
      
      ntox.st = sum(tox)
      
      t.tox[tox == 0] = Tobs
      
    }
    return(list(
      tox = tox,
      t.tox = t.tox,
      ntox.st = ntox.st
    ))
    
  }
  
  
  dselect = rep(0, ntrial)
  ndose = length(teff)
  N = matrix(rep(0, ndose * ntrial), ncol = ndose)
  YTOX = matrix(rep(0, ndose * ntrial), ncol = ndose)
  YEFF = matrix(rep(0, ndose * ntrial), ncol = ndose)
  durationV = rep(0, ntrial)
  
  
  for (trial in 1:ntrial) {
    yeff.i = NULL  #efficacy indicator for each subject
    ytox.i = NULL
    
    dv = NULL
    #dose for each subject
    
    ytox = rep(0, ndose) # number of toxicity at each dose by interim time
    yeff = rep(0, ndose) # number of efficacy at each dose by interim time
    
    n.eff = rep(0, ndose)    # Ess of efficacy at each dose
    n.tox = rep(0, ndose)    # Ess of toxicity at each dose
    
    ytox.d = rep(0, ndose)
    yeff.d = rep(0, ndose)
    n.d = rep(0, ndose)
    
    efft.enter = NULL
    efft.event = NULL
    efft.decision = 0
    
    toxt.enter = NULL
    toxt.event = NULL
    toxt.decision = 0
    
    d = 1 ##dose start from level 1
    
    for (i in 1:ncohort) {
      if(d==0){break}
      # generate data for the new patient
      for (j in 1:cohortsize) {
        if (j == 1) {
          efft.enter = c(efft.enter, efft.decision)
          toxt.enter = efft.enter
        }
        else {
          if (dist3 == 1) {
            efft.enter = c(efft.enter, efft.enter[length(efft.enter)] + runif(1, 0, 2 /
                                                                                accrual))
            toxt.enter = efft.enter
          }
          if (dist3 == 2) {
            efft.enter = c(efft.enter, efft.enter[length(efft.enter)] + rexp(1, rate =
                                                                               accrual))
            toxt.enter = efft.enter
          }
        }
      }
      
      obscohort.eff = gen.tite(dist1,
                               cohortsize,
                               teff[d],
                               alpha = alpha1,
                               Tobs =
                                 maxt1)
      obscohort.tox = gen.tite(dist2,
                               cohortsize,
                               ttox[d],
                               alpha = alpha2,
                               Tobs =
                                 maxt2)
      
      efft.event = c(efft.event, obscohort.eff$t.tox)
      toxt.event = c(toxt.event, obscohort.tox$t.tox)
      
      yeff.i = c(yeff.i, obscohort.eff$tox)
      ytox.i = c(ytox.i, obscohort.tox$tox)
      
      dv = c(dv, rep(d, cohortsize))
      
      efft.decision = efft.enter[length(efft.enter)]
      toxt.decision = toxt.enter[length(toxt.enter)]
      
      pending = 1
      
      while (pending == 1) {
        pending = 0
        
        if (i == ncohort) {
          efft.decision = efft.decision + maxt1
          toxt.decision = toxt.decision + maxt2
        }
        else {
          if (dist3 == 1) {
            efft.decision = efft.decision + runif(1, 0, 2 / accrual)
            toxt.decision = efft.decision 
          }
          if (dist3 == 2) {
            efft.decision = efft.decision + rexp(1, rate = accrual)
            toxt.decision = efft.decision
          }
        }
        
        # determine which observation are observed
        delta.eff = ((efft.enter + efft.event) <= efft.decision)
        delta.tox = ((toxt.enter + toxt.event) <= toxt.decision)
        
        t.eff = pmin(efft.event, efft.decision - efft.enter, maxt1)
        ## used for recording potential censoring time
        t.tox = pmin(toxt.event, toxt.decision - toxt.enter, maxt2)
        ## used for recording potential censoring time
        
        ######define whole data
        for (dd in 1:ndose) {
          cset1 = (dv == dd)
          delta.eff.curr1 = delta.eff[cset1]
          delta.tox.curr1 = delta.tox[cset1]
          
          t.eff.curr1 = t.eff[cset1]
          t.tox.curr1 = t.tox[cset1]
          
          neff.curr1 = sum((yeff.i[cset1])[delta.eff.curr1 == 1])
          ntox.curr1 = sum((ytox.i[cset1])[delta.tox.curr1 == 1])
          
          totalt.eff1 = t.eff.curr1[delta.eff.curr1 == 0]
          totalt.tox1 = t.tox.curr1[delta.tox.curr1 == 0]
          
          totalt.eff1 = sum(totalt.eff1) / maxt1
          totalt.tox1 = sum(totalt.tox1) / maxt2
          
          n.curr1 = sum(cset1)
          
          n.eff.pend1 = sum(delta.eff[cset1] == 0)
          n.tox.pend1 = sum(delta.tox[cset1] == 0)
          
          effobs1 = sum(delta.eff[cset1])
          toxobs1 = sum(delta.tox[cset1])
          
          
          n.eff[dd] = n.curr1 - n.eff.pend1 + totalt.eff1
          yeff[dd] = neff.curr1
          
          n.tox[dd] = n.curr1 - n.tox.pend1 + totalt.tox1
          ytox[dd] = ntox.curr1
        }
        ####current dose data
        cset = (dv == d)
        
        delta.eff.curr = delta.eff[cset]
        delta.tox.curr = delta.tox[cset]
        
        t.eff.curr = t.eff[cset]
        t.tox.curr = t.tox[cset]
        
        neff.curr = sum((yeff.i[cset])[delta.eff.curr == 1])
        ntox.curr = sum((ytox.i[cset])[delta.tox.curr == 1])
        
        totalt.eff = t.eff.curr[delta.eff.curr == 0]
        totalt.tox = t.tox.curr[delta.tox.curr == 0]
        
        totalt.eff = sum(totalt.eff) / maxt1
        totalt.tox = sum(totalt.tox) / maxt2
        
        n.curr = sum(cset)
        
        n.eff.pend = sum(delta.eff[cset] == 0)
        n.tox.pend = sum(delta.tox[cset] == 0)
        
        effobs = sum(delta.eff[cset])
        toxobs = sum(delta.tox[cset])
        
        
        n.eff[d] = n.curr - n.eff.pend + totalt.eff
        yeff[d] = neff.curr
        
        n.tox[d] = n.curr - n.tox.pend + totalt.tox
        ytox[d] = ntox.curr
        
        if ((effobs > (0.5*n.curr)) & (toxobs > (0.5*n.curr))){
          pending=0
          AT = adm_tox(n = n.tox, ytox = ytox)
          AR = adm_eff(n = n.eff, yeff = yeff)
          try = length(n.eff[n.eff != 0])
          if ((try %in% AT) & (try < ndose)) {
            d = d + 1
          } else{
            A = intersect(AT, AR)
            if (length(A) == 0) {
              d = 0
              break
            } else{
              OBD = A[which.min(AIC(n = n.eff, yeff = yeff)[A])]
              if (OBD > d) {
                d = d + 1
              } else if (OBD < d) {
                d = d - 1
              } else {
                d = d
              }
            }
          }
        }
        else{
          pending = 1
        }
      }
    }
    
    for (k in 1:ndose) {
      ytox.d[k] = sum(ytox.i[dv == k])
      yeff.d[k] = sum(yeff.i[dv == k])
      n.d[k] = sum(dv == k)
    }
    
    if (d == 0) {
      dselect[trial] = 0
      N[trial, ] = n.d
      YTOX[trial, ] = ytox.d
      YEFF[trial, ] = yeff.d
      durationV[trial] = max(toxt.decision,efft.decision)
    } else{
      AT = adm_tox(n = n.tox, ytox = ytox)
      AR = adm_eff(n = n.eff, yeff = yeff)
      A = intersect(AT, AR)
      if (length(A) == 0){
        dselect[trial] = 0
        N[trial, ] = n.d
        YTOX[trial, ] = ytox.d
        YEFF[trial, ] = yeff.d
        durationV[trial] = max(toxt.decision,efft.decision)
      }
      else{
        dselect[trial] = A[which.min(AIC(n = n.eff, yeff = yeff)[A])]
        N[trial, ] = n.d
        YTOX[trial, ] = ytox.d
        YEFF[trial, ] = yeff.d
        durationV[trial] = max(toxt.decision,efft.decision)
      }
    }
  }
  
  selpercent = rep(0, ndose + 1)
  patpercent = matrix(rep(0, ntrial * ndose), ncol = ntrial, nrow = ndose)
  efficacy = rep(0, ntrial)
  toxicity = rep(0, ntrial)
  f <- function(x) {
    x[i] / sum(x)
  }
  
  ## Summarize results
  for (i in 0:ndose) {
    selpercent[(i + 1)] = sum(dselect == i) / ntrial * 100
  }
  print("selection probablity")
  cat(formatC(selpercent, digits = 1, format = "f"), sep = " ", "\n")
  for (i in 1:ndose) {
    patpercent[i, ] = apply(N, 1, f)
  }
  print("average percent of patients")
  cat(formatC(
    apply(patpercent, 1, mean) * 100,
    digits = 1,
    format = "f"
  ),
  sep = " ", "\n")
  print("average number of patients")
  cat(formatC(c(apply(N, 2, mean), sum(
    apply(N, 2, mean)
  )), digits = 1, format = "f"),
  sep = " ", "\n")
  print("average number of patients response to efficacy")
  cat(formatC(c(apply(YEFF, 2, mean), sum(
    apply(YEFF, 2, mean)
  )), digits = 1, format = "f"),
  sep = " ", "\n")
  for (i in 1:ntrial) {
    efficacy[i] = sum(YEFF[i, ]) / sum(N[i, ])
  }
  print("average percent of efficacy")
  cat(formatC(mean(efficacy) * 100, digits = 1, format = "f"),
      sep = " ", "\n")
  print("average number of patients response to toxicity")
  cat(formatC(c(apply(YTOX, 2, mean), sum(
    apply(YTOX, 2, mean)
  )), digits = 1, format = "f"),
  sep = " ", "\n")
  for (i in 1:ntrial) {
    toxicity[i] = sum(YTOX[i, ]) / sum(N[i, ])
  }
  print("average percent of toxicity")
  cat(formatC(mean(toxicity) * 100, digits = 1, format = "f"),
      sep = " ", "\n")
  print("average trial duration")
  cat(formatC(mean(durationV) , digits = 1, format = "f"),
      sep = " ", "\n")
}


############################################################################
## Function for real trial dose-finding procedure of mISO-B
##
## ndose: number of dose levels 
## eff_event: efficacy event for each individual by interim time
## tox_event: toxicity event for each individual by interim time
## eff_time: time from entering the trial to efficacy event, 
##           if do not experience the efficacy by interim time, simply plug in the time assessment window for efficacy
## tox_time: time from entering the trial to toxicity event, 
##           if do not experience the toxicity by interim time, simply plug in the time assessment window for toxicity
## enter_time: time for individuals entering the trial
## decision_time: interim decison time
## d: dose level for each individual
## maxt1: time assess window for efficacy
## maxt2: time assess window for toxicity
############################################################################

mISO_B.df = function(ndose,eff_event,tox_event,eff_time,tox_time,enter_time,decision_time,
                     d,phiT=0.3,muT=0.9,phiE=0.5,muE=0.85,
                     alphaT=0.5,betaT=0.5,alphaE=0.5,betaE=0.5,maxt1,maxt2){
  library("Iso")
  
  ###find admissible set of toxicity
  adm_tox <- function(n, ytox) {
    nn = n[n != 0]
    yytox = ytox[which(n != 0)]
    at = alphaT + yytox
    bt = betaT + nn - yytox
    Tox_prob = 1 - pbeta(phiT, at, bt)
    AT_naive = which(Tox_prob < muT)
    if (length(AT_naive)==0){
      AT=AT_naive
    }
    else{
      full_seq = seq(min(AT_naive),max(AT_naive),1)
      if (length(setdiff(full_seq, AT_naive)) == 0) {
        AT = AT_naive
      }
      else{
        AT = AT_naive[AT_naive < min(setdiff(full_seq, AT_naive))]
      }
    }
    return(AT)
  }
  
  ###find admissible set of efficacy
  adm_eff <- function(n, yeff) {
    nn = n[n != 0]
    yyeff = yeff[which(n != 0)]
    ae = alphaE + yyeff
    be = betaE + nn - yyeff
    Eff_prob = pbeta(phiE, ae, be)
    AR_naive = which(Eff_prob < muE)
    if (length(AR_naive)==0){
      AR=AR_naive
    }
    else{
      full_seq = seq(min(AR_naive),max(AR_naive),1)
      if (length(setdiff(full_seq, AR_naive)) == 0) {
        AR = AR_naive
      } else {
        AR = AR_naive[AR_naive > max(setdiff(full_seq, AR_naive))]
      }
    }
    return(AR)
  }
  
  
  ###AIC model selection
  
  AIC <- function(n, yeff) {
    aic = rep(0, length(n[n != 0]))
    for (l in 1:length(n[n != 0])) {
      if (l == 1) {
        jstar = length(n[n != 0])
        sm = sum(yeff[1:jstar])
        nm = sum(n[1:jstar])
        ql = sm / nm
        qfit = pava(y = ql, w = nm)
        ql_hat = qfit[1]
        likhood = (ql_hat ^ sm) * (1 - ql_hat) ^ (nm - sm)
        aic[l] = 2 * length(qfit) - 2 * log(likhood)
      } else{
        jstar = length(n[n != 0])
        sk = (yeff[which(n != 0)])[1:(l - 1)]
        nk = (n[n != 0])[1:(l - 1)]
        qk = sk / nk
        sm = sum(yeff[l:jstar])
        nm = sum(n[l:jstar])
        ql = sm / nm
        qfit = pava(y = c(qk, ql), w = c(nk, nm))
        qk_hat = qfit[1:(l - 1)]
        ql_hat = qfit[l]
        likhood = prod((qk_hat ^ sk) * (1 - qk_hat) ^ (nk - sk)) * (ql_hat ^
                                                                      sm) * (1 - ql_hat) ^ (nm - sm)
        aic[l] = 2 * length(qfit) - 2 * log(likhood)
      }
    }
    return(aic)
  }
  
  delta.eff<-((enter_time + eff_time) <= (decision_time-1))
  delta.tox<-((enter_time + tox_time) <= (decision_time-1))
  eff_t<-pmin(eff_time,(decision_time-1)-enter_time,90)
  tox_t<-pmin(tox_time,(decision_time-1)-enter_time,90)
  dv = d

  ytox = rep(0, ndose) # number of toxicity at each dose by interim time
  yeff = rep(0, ndose) # number of efficacy at each dose by interim time
  n.eff = rep(0, ndose)    # Ess of efficacy at each dose
  n.tox = rep(0, ndose)    # Ess of toxicity at each dose
  
  for (dd in 1:ndose) {
    cset1 = (dv == dd)
    delta.eff.curr1 = delta.eff[cset1]
    delta.tox.curr1 = delta.tox[cset1]
    
    t.eff.curr1 = eff_t[cset1]
    t.tox.curr1 = tox_t[cset1]
         
    neff.curr1 = sum((eff_event[cset1])[delta.eff.curr1 == 1])
    ntox.curr1 = sum((tox_event[cset1])[delta.tox.curr1 == 1])
             
    totalt.eff1 = t.eff.curr1[delta.eff.curr1 == 0]
    totalt.tox1 = t.tox.curr1[delta.tox.curr1 == 0]
              
    totalt.eff1 = sum(totalt.eff1) / maxt1
    totalt.tox1 = sum(totalt.tox1) / maxt2
            
    n.curr1 = sum(cset1)
                   
    n.eff.pend1 = sum(delta.eff[cset1] == 0)
    n.tox.pend1 = sum(delta.tox[cset1] == 0)
                    
    effobs1 = sum(delta.eff[cset1])
    toxobs1 = sum(delta.tox[cset1])
                       
    n.eff[dd] = n.curr1 - n.eff.pend1 + totalt.eff1
    yeff[dd] = neff.curr1
    n.tox[dd] = n.curr1 - n.tox.pend1 + totalt.tox1
    ytox[dd] = ntox.curr1
  }
  
  ####current dose data
  cset = (dv == tail(d,1))
  
  delta.eff.curr = delta.eff[cset]
  delta.tox.curr = delta.tox[cset]
  
  t.eff.curr = eff_t[cset]
  t.tox.curr = tox_t[cset]
  
  neff.curr = sum((eff_event[cset])[delta.eff.curr == 1])
  ntox.curr = sum((tox_event[cset])[delta.tox.curr == 1])
  
  totalt.eff = t.eff.curr[delta.eff.curr == 0]
  totalt.tox = t.tox.curr[delta.tox.curr == 0]
  
  totalt.eff = sum(totalt.eff) / maxt1
  totalt.tox = sum(totalt.tox) / maxt2
  
  n.curr = sum(cset)
  
  effobs = sum(delta.eff[cset])
  toxobs = sum(delta.tox[cset])
  
  
  AT = adm_tox(n = n.tox, ytox = ytox)
  AR = adm_eff(n = n.eff, yeff = yeff)
  try = length(n.eff[n.eff != 0])
  
  nd=tail(d,1)
  
  if ((effobs > (0.5*n.curr)) & (toxobs > (0.5*n.curr))){
    AT = adm_tox(n = n.tox, ytox = ytox)
    AR = adm_eff(n = n.eff, yeff = yeff)
    try = length(n.eff[n.eff != 0])
    if ((try %in% AT) & (try < ndose)) {
      nd = nd + 1
    } else{
      A = intersect(AT, AR)
      if (length(A) == 0) {
        nd = 0
      } else{
        OBD = A[which.min(AIC(n = n.eff, yeff = yeff)[A])]
        if (OBD > nd) {
          nd = nd + 1
        } else if (OBD < nd) {
          nd = nd - 1
        } else {
          nd = nd
        }
      }
    }
    return(list("next dose"=nd)) 
  }
  else{
    return(list("next dose"="suspend the trial")) 
  }

}

########################################## an example #######################################

####### generate operating characteristics of the designs ##########
# assume that the target toxicity rate is 30%, the target efficacy rate is 50%, 
# and the sample size is 60 in cohorts of size 3
eff=c(0.4,0.6,0.6,0.6,0.6,0.6)
tox=c(0.03,0.1,0.2,0.3,0.4,0.5)
mISO.oc(
  ttox = tox,
  teff = eff,
  cohortsize = 3,
  ncohort = 20,
  phiT = 0.3,
  muT = 0.9,
  phiE = 0.5,
  muE = 0.85,
  alphaT = 0.5,
  betaT = 0.5,
  alphaE = 0.5,
  betaE = 0.5,
  ntrial = 10000
)

#[1] "selection probablity"
#14.9 14.7 54.1 10.2 4.4 1.5 0.1 
#[1] "average percent of patients"
#25.5 36.6 15.1 11.4 8.1 3.3 
#[1] "average number of patients"
#13.7 20.8 8.3 6.1 4.2 1.6 54.6 
#[1] "average number of patients response to efficacy"
#5.5 12.5 5.0 3.6 2.5 0.9 30.0 
#[1] "average percent of efficacy"
#53.8 
#[1] "average number of patients response to toxicity"
#0.4 2.1 1.6 1.8 1.7 0.8 8.4 
#[1] "average percent of toxicity"
#16.4 


# Assume that the target toxicity rate is 30%, the target efficacy rate is 50%, and the sample size is 60 in cohorts of size 3. 
# Assessment windows are both 3 months, and patients are enrolled at the rate of 3 patients per month. 
# The underlying distribution to efficacy and toxicity is weibull distribution by controlling that 
#    50% of events occur in the latter half of the assessment window and the patient accrual time is uniformly distributed. 

mISO_B.oc(
  ttox = tox,
  teff = eff,
  cohortsize = 3,
  ncohort = 20,
  phiT = 0.3,
  muT = 0.9,
  phiE = 0.5,
  muE = 0.85,
  alphaT = 0.5,
  betaT = 0.5,
  alphaE = 0.5,
  betaE = 0.5,
  ntrial = 10000,
  accrual = 3,
  maxt1 = 3,
  maxt2 = 3,
  alpha1 = 0.5,
  alpha2 = 0.5,
  dist1 = 2,
  dist2 = 2,
  dist3 = 1
)

#[1] "selection probablity"
#14.9 16.4 53.4 9.6 4.2 1.3 0.1 
#[1] "average percent of patients"
#26.8 34.5 15.3 11.6 8.2 3.5 
#[1] "average number of patients"
#14.5 19.6 8.4 6.2 4.3 1.8 54.8 
#[1] "average number of patients response to efficacy"
#5.8 11.8 5.0 3.7 2.6 1.0 29.9 
#[1] "average percent of efficacy"
#53.6 
#[1] "average number of patients response to toxicity"
#0.4 2.0 1.7 1.9 1.7 0.9 8.5 
#[1] "average percent of toxicity"
#16.6 
#[1] "average trial duration"
#39.7 

#############  trial conduct  #############
# mISO
# assume that the first cohort is treated at dose level 1 without toxicity and one efficacy response, yielding the data
n=c(3,0,0,0,0)
ytox=c(0,0,0,0,0)
yeff=c(1,0,0,0,0)
d=1
## call dose-finding function to obtain the recommended dose assignment for the next cohort
mISO.df(n=n,ytox = ytox,yeff=yeff,d=d)
##$`next dose`
##[1] 2


## Hence, we escalate the dose level and treat the second cohort at dose level 2
## assume that in cohort 2, no patient experience toxicity and 2 patients have efficacy response at dose level 2, yielding the updated data
n=c(3,3,0,0,0)
ytox=c(0,0,0,0,0)
yeff=c(1,2,0,0,0)
d=2
## call dose-finding function to obtain the recommmended dose assignment for the next cohort
mISO.df(n=n,ytox = ytox,yeff=yeff,d=d)
##$`next dose`
##[1] 3


## Hence, we escalate the dose level and treat the thrid cohort at dose level 3
.
.
## ....... repeat the above procedure for incoming cohorts................
.
.
## Suppose after the last cohort of patients have been treated, the updated data are
n=c(3,6,12,6,3)
ytox=c(0,1,2,2,2)
yeff=c(0,2,5,3,2)

## call dose-recommendation function to obtain the recommended optimal biological dose (OBD)
mISO.select(n=n,ytox = ytox,yeff = yeff)
##$OBD
##[1] 2


# mISO-B
## using the example in trial implementation 
## call dose-finding function to obtain the recommended dose assignment for the next cohort at day 102
mISO_B.df(ndose=5,eff_event = c(0,0,0),tox_event = c(0,0,0),eff_time = c(90,90,90),tox_time = c(90,90,90),
          enter_time = c(1,11,21),decision_time = c(102,102,102),d=c(1,1,1),maxt1 = 90,maxt2 = 90)
#$`next dose`
#[1] 2

## call dose-finding function to obtain the recommended dose assignment for the next cohort at day 203
mISO_B.df(ndose=5,eff_event = c(0,0,0,0,1,0),tox_event = c(0,0,0,0,0,0),eff_time = c(90,90,90,90,50,90),
          tox_time = c(90,90,90,90,90,90),enter_time = c(1,11,21,102,112,122),decision_time = c(203,203,203,203,203,203),
          d=c(1,1,1,2,2,2),maxt1 = 90,maxt2 = 90)
#$`next dose`
#[1] 3

## suppose we want to obtain the recommended dose assignment for the next cohort at day 193,
## we need to supsend the trail until 50% of the patients in the current dose do not have pending data 
mISO_B.df(ndose=5,eff_event = c(0,0,0,0,1,0),tox_event = c(0,0,0,0,0,0),eff_time = c(90,90,90,90,50,90),
          tox_time = c(90,90,90,90,90,90),enter_time = c(1,11,21,102,112,122),decision_time = c(193,193,193,193,193,193),
          d=c(1,1,1,2,2,2),maxt1 = 90,maxt2 = 90)

## ....... repeat the above procedure for incoming cohorts................


#$`next dose`
#[1] "suspend the trial"

## Suppose after the last cohort of patients have been treated, the updated data are
n=c(3,3,3,6,3)
ytox = c(0,0,1,5,2)
yeff = c(0,1,1,4,2)

## call dose-recommendation function to obtain the recommended optimal biological dose (OBD)
mISO.select(n=n,ytox = ytox,yeff = yeff)
##$OBD
##[1] 2
