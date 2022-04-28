mISO_B <- function(ttox,
                     teff,
                     cohortsize = 3,
                     ncohort = 20,
                     ntrial = 5000,
                     thetaT = 0.3,
                     lambdaT = 0.95,
                     thetaR = 0.5,
                     lambdaR = 0.8,
                     dist1 = 1,# distribution of efficacy
                     dist2 = 1,#distribution of toxicity
                     dist3 = 1,#patient arrival
                     alpha = 0.5,
                     maxt1 = 3,# time assess window for efficacy
                     maxt2=3,# time assess window for toxicity
                     accrual = 2) {
  library("Iso")
  
  ###find admissible set of toxicity
  adm_tox <- function(n, ytox) {
    nn = n[n != 0]
    yytox = ytox[which(n != 0)]
    at = 0.5 + yytox
    bt = 0.5 + nn - yytox
    Tox_prob = 1 - pbeta(thetaT, at, bt)
    AT_naive = which(Tox_prob < lambdaT)
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
    ae = 0.5 + yyeff
    be = 0.5 + nn - yyeff
    Eff_prob = pbeta(thetaR, ae, be)
    AR_naive = which(Eff_prob < lambdaR)
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
                               alpha = alpha,
                               Tobs =
                                 maxt1)
      obscohort.tox = gen.tite(dist2,
                               cohortsize,
                               ttox[d],
                               alpha = alpha,
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
####example
eff = c(0.4, 0.6, 0.6, 0.6, 0.6, 0.6)
tox = c(0.03, 0.1, 0.2, 0.3, 0.4, 0.5)
mISO_B(
  ttox = tox,
  teff = eff,
  cohortsize = 3,
  ncohort = 20,
  thetaT = 0.3,
  lambdaT = 0.9,
  thetaR = 0.5,
  lambdaR = 0.85,
  ntrial = 10000,
  accrual = 3,
  maxt1 = 3,
  maxt2 = 3,
  alpha = 0.5,
  dist1 = 2,
  dist2 = 2,
  dist3 = 1
  
)
