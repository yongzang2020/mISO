#####
mISO = function(ttox,
                       teff,
                       cohortsize = 3,
                       ncohort = 10,
                       ntrial = 5000,
                       thetaT = 0.3,
                       lambdaT = 0.95,
                       thetaR = 0.5,
                       lambdaR = 0.7) {
  library("Iso")
  
  ###find admissible set of toxicity
  adm_tox <- function(n, ytox) {
    nn = n[n != 0]
    yytox = ytox[which(n != 0)]
    at = 0.5 + yytox
    bt = 0.5 + nn - yytox
    Tox_prob = 1 - pbeta(thetaT, at, bt)
    AT_naive = which(Tox_prob < lambdaT)
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
    ae = 0.5 + yyeff
    be = 0.5 + nn - yyeff
    Eff_prob = pbeta(thetaR, ae, be)
    AR_naive = which(Eff_prob < lambdaR)
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

#######example
eff=c(0.4,0.6,0.6,0.6,0.6,0.6)
tox=c(0.03,0.1,0.2,0.3,0.4,0.5)
mISO(
  ttox = tox,
  teff = eff,
  cohortsize = 3,
  ncohort = 20,
  thetaT = 0.3,
  lambdaT = 0.9,
  thetaR = 0.5,
  lambdaR = 0.85,
  ntrial = 10000
)
