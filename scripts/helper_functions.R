
fun_generate_priors <- function(analyte_type, logP=2.2){
  # Priors based on: Comparison of Chromatographic Stationary Phases Using a Bayesian-Based Multilevel Model
  # Paweł Wiczling and Agnieszka Kamedulska
  # Analytical Chemistry 2024 96 (3), 1310-1319
  # DOI: 10.1021/acs.analchem.3c04697
  
  param <- list(
  logkwHat = 3.60,  
  S1mHat  = 4.92,   
  dS1Hat  = 0.61,   
  dlogkwHatA = -0.79, 
  dlogkwHatB = -0.97,  
  dS1mHatA = 0.17,   
  dS1mHatB = 0.12,   
  ddS1HatA = 0.28,   
  ddS1HatB = -0.67, 
  logS2mHat = -0.31, 
  dlogS2Hat =  0.42,   
  beta_logkw = 0.84,   
  beta_S1= 0.51,  
  dlogkTHat = -0.09, 
  apHA= -0.03,  
  apHB= 0.08,  
  omega_logkw = 0.92,   
  omega_S1m = 0.93,   
  omega_dSa = 0.55,   
  omegaT = 0.03,  
  kappa_dlogkwHat = 0.59 ,  
  kappa_dS1m = 0.69,   
  kappa_ddS1 = 0.55,  
  rho_logkw_S1m = 0.87,  
  msigma = 0.39,  
  ssigma = 0.81)

  
  if (analyte_type == "N_MeOH") {
    p=2 
    miu = c(param$logkwHat, param$S1mHat)
    beta  = c(param$beta_logkw,param$beta_S1)
    theta = miu + beta * (logP-2.2)
    msigma_sqaure = param$msigma^2
    omega = c(param$omega_logkw, param$omega_S1m)
    rho = matrix(c(1,param$rho_logkw_S1m,
                   param$rho_logkw_S1m,1),2,2,byrow = TRUE)
    S = diag(omega)%*% rho %*%diag(omega)
    S_inv = inv(S)
  }
  
  
  if (analyte_type == "N_ACN") {
    p=2 
    miu = c(param$logkwHat, param$S1mHat+param$dS1Hat)
    beta  = c(param$beta_logkw,param$beta_S1)
    theta = miu + beta * (logP-2.2)
    msigma_sqaure = param$msigma^2
    omega = c(param$omega_logkw, param$omega_S1m, param$omega_dSa)
    rho = matrix(c(1,param$rho_logkw_S1m,0,
                   param$rho_logkw_S1m,1,0,
                   0,0,1),3,3,byrow = TRUE)
    tS = diag(omega)%*% rho %*%diag(omega)
    
    B = matrix(c(1,0,0,
                 0,1,1),2,3,byrow = TRUE)
    S= B%*%tS%*%t(B)
    S_inv = inv(S)
  }
  
  if (analyte_type == "N") {
  p=3 
  miu = c(param$logkwHat, param$S1mHat, param$dS1Hat)
  beta  = c(param$beta_logkw,param$beta_S1,0)
  theta = miu + beta * (logP-2.2)
  msigma_sqaure = param$msigma^2
  omega = c(param$omega_logkw, param$omega_S1m, param$omega_dSa)
  rho = matrix(c(1,param$rho_logkw_S1m,0,
                   param$rho_logkw_S1m,1,0,
                   0,0,1),3,3,byrow = TRUE)
  S = diag(omega)%*% rho %*%diag(omega)
  S_inv = inv(S)
    
  }
  
  if (analyte_type == "A") {
    p=3
    beta  = c(param$beta_logkw,param$beta_S1,0,0,0,0)
    miu = c(param$logkwHat, param$S1mHat, param$dS1Hat, param$dlogkwHatA, param$dS1mHatA, param$ddS1HatA)+beta*(logP-2.2)
    msigma_sqaure = param$msigma^2
    omega = c(param$omega_logkw, param$omega_S1m, param$omega_dSa,param$kappa_dlogkwHat,param$kappa_dS1m,param$kappa_ddS1)
    rho = diag(rep(1,6))
    rho[1,2] = param$rho_logkw_S1m
    rho[2,1] = param$rho_logkw_S1m  
    tS = diag(omega)%*% rho %*%diag(omega)
    
    B = matrix(c(1,0,0,1,0,0,
                0,1,0,0,1,0,
                0,0,1,0,0,1),3,6,byrow = TRUE)

    theta= B%*%miu
    S = B%*%tS%*%t(B) 
    S_inv = inv(S)
 
  }
  
  if (analyte_type == "A_MeOH") {
    p=2
    beta  = c(param$beta_logkw,param$beta_S1,0,0,0,0)
    miu = c(param$logkwHat, param$S1mHat, param$dS1Hat, param$dlogkwHatA, param$dS1mHatA, param$ddS1HatA)+beta*(logP-2.2)
    msigma_sqaure = param$msigma^2
    omega = c(param$omega_logkw, param$omega_S1m, param$omega_dSa,param$kappa_dlogkwHat,param$kappa_dS1m,param$kappa_ddS1)
    rho = diag(rep(1,6))
    rho[1,2] = param$rho_logkw_S1m
    rho[2,1] = param$rho_logkw_S1m  
    tS = diag(omega)%*% rho %*%diag(omega)
    
    B = matrix(c(1,0,0,1,0,0,
                 0,1,0,0,1,0),2,6,byrow = TRUE)
    
    theta= B%*%miu
    S = B%*%tS%*%t(B) 
    S_inv = inv(S)
    
  }
  
  
  if (analyte_type == "A_ACN") {
    p=2
    beta  = c(param$beta_logkw,param$beta_S1,0,0,0,0)
    miu = c(param$logkwHat, param$S1mHat, param$dS1Hat, param$dlogkwHatA, param$dS1mHatA, param$ddS1HatA)+beta*(logP-2.2)
    msigma_sqaure = param$msigma^2
    omega = c(param$omega_logkw, param$omega_S1m, param$omega_dSa,param$kappa_dlogkwHat,param$kappa_dS1m,param$kappa_ddS1)
    rho = diag(rep(1,6))
    rho[1,2] = param$rho_logkw_S1m
    rho[2,1] = param$rho_logkw_S1m  
    tS = diag(omega)%*% rho %*%diag(omega)
    
    B = matrix(c(1,0,0,1,0,0,
                 0,1,1,0,1,1),2,6,byrow = TRUE)
    
    theta= B%*%miu
    S = B%*%tS%*%t(B) 
    S_inv = inv(S)
    
  }
  
  if (analyte_type == "NA") {
    p=6
    beta  = c(param$beta_logkw,param$beta_S1,0,0,0,0)
    miu = c(param$logkwHat, param$S1mHat, param$dS1Hat, param$dlogkwHatA, param$dS1mHatA, param$ddS1HatA)+beta*(logP-2.2)
    msigma_sqaure = param$msigma^2
    omega = c(param$omega_logkw, param$omega_S1m, param$omega_dSa,param$kappa_dlogkwHat,param$kappa_dS1m,param$kappa_ddS1)
    rho = diag(rep(1,6))
    rho[1,2] = param$rho_logkw_S1m
    rho[2,1] = param$rho_logkw_S1m  
    tS = diag(omega)%*% rho %*%diag(omega)
    theta = miu
    S = tS
    S_inv = inv(S)

  }
  
  if (analyte_type == "B") {
    p=3
    beta  = c(param$beta_logkw,param$beta_S1,0,0,0,0)
    miu = c(param$logkwHat, param$S1mHat, param$dS1Hat, param$dlogkwHatB, param$dS1mHatB, param$ddS1HatB)+beta*(logP-2.2)
    msigma_sqaure = param$msigma^2
    omega = c(param$omega_logkw, param$omega_S1m, param$omega_dSa,param$kappa_dlogkwHat,param$kappa_dS1m,param$kappa_ddS1)
    rho = diag(rep(1,6))
    rho[1,2] = param$rho_logkw_S1m
    rho[2,1] = param$rho_logkw_S1m  
    tS = diag(omega)%*% rho %*%diag(omega)
    
    B = matrix(c(1,0,0,1,0,0,
                 0,1,0,0,1,0,
                 0,0,1,0,0,1),3,6,byrow = TRUE)
    
    theta= B%*%miu
    S = B%*%tS%*%t(B) 
    S_inv = inv(S)
  }  
  
  if (analyte_type == "B_MeOH") {
    p=2
    beta  = c(param$beta_logkw,param$beta_S1,0,0,0,0)
    miu = c(param$logkwHat, param$S1mHat, param$dS1Hat, param$dlogkwHatB, param$dS1mHatB, param$ddS1HatB)+beta*(logP-2.2)
    msigma_sqaure = param$msigma^2
    omega = c(param$omega_logkw, param$omega_S1m, param$omega_dSa,param$kappa_dlogkwHat,param$kappa_dS1m,param$kappa_ddS1)
    rho = diag(rep(1,6))
    rho[1,2] = param$rho_logkw_S1m
    rho[2,1] = param$rho_logkw_S1m  
    tS = diag(omega)%*% rho %*%diag(omega)
    
    B = matrix(c(1,0,0,1,0,0,
                 0,1,0,0,1,0),2,6,byrow = TRUE)
    
    theta= B%*%miu
    S = B%*%tS%*%t(B) 
    S_inv = inv(S)
  } 
  
  if (analyte_type == "B_ACN") {
    p=2
    beta  = c(param$beta_logkw,param$beta_S1,0,0,0,0)
    miu = c(param$logkwHat, param$S1mHat, param$dS1Hat, param$dlogkwHatB, param$dS1mHatB, param$ddS1HatB)+beta*(logP-2.2)
    msigma_sqaure = param$msigma^2
    omega = c(param$omega_logkw, param$omega_S1m, param$omega_dSa,param$kappa_dlogkwHat,param$kappa_dS1m,param$kappa_ddS1)
    rho = diag(rep(1,6))
    rho[1,2] = param$rho_logkw_S1m
    rho[2,1] = param$rho_logkw_S1m  
    tS = diag(omega)%*% rho %*%diag(omega)
    
    B = matrix(c(1,0,0,1,0,0,
                 0,1,1,0,1,1),2,6,byrow = TRUE)
    
    theta= B%*%miu
    S = B%*%tS%*%t(B) 
    S_inv = inv(S)
  } 
  
  if (analyte_type == "NB") {
    p=6
    beta  = c(param$beta_logkw,param$beta_S1,0,0,0,0)
    miu = c(param$logkwHat, param$S1mHat, param$dS1Hat, param$dlogkwHatB, param$dS1mHatB, param$ddS1HatB)+beta*(logP-2.2)
    msigma_sqaure = param$msigma^2
    omega = c(param$omega_logkw, param$omega_S1m, param$omega_dSa,param$kappa_dlogkwHat,param$kappa_dS1m,param$kappa_ddS1)
    rho = diag(rep(1,6))
    rho[1,2] = param$rho_logkw_S1m
    rho[2,1] = param$rho_logkw_S1m  
    tS = diag(omega)%*% rho %*%diag(omega)
    theta = miu
    S = tS  
    S_inv = inv(S)
  }  
  
  priors = list(
    p=p,
    beta=beta,
    theta = theta,
    msigma_sqaure = msigma_sqaure,
    S = S,
    S_inv = S_inv)
  
  return(priors)
}

fun_generate_conditional_priors <- function(analyte_type, logP=2.2){
  # Priors based on: Comparison of Chromatographic Stationary Phases Using a Bayesian-Based Multilevel Model
  # Paweł Wiczling and Agnieszka Kamedulska
  # Analytical Chemistry 2024 96 (3), 1310-1319
  # DOI: 10.1021/acs.analchem.3c04697
  # logk = 1 for fix = 0.4 in MeOH
  param <- list(
    logkwHat = 3.60,  
    S1mHat  = 4.92,   
    dS1Hat  = 0.61,   
    dlogkwHatA = -0.79, 
    dlogkwHatB = -0.97,  
    dS1mHatA = 0.17,   
    dS1mHatB = 0.12,   
    ddS1HatA = 0.28,   
    ddS1HatB = -0.67, 
    logS2mHat = -0.31, 
    dlogS2Hat =  0.42,   
    beta_logkw = 0.84,   
    beta_S1= 0.51,  
    dlogkTHat = -0.09, 
    apHA= -0.03,  
    apHB= 0.08,  
    omega_logkw = 0.92,   
    omega_S1m = 0.93,   
    omega_dSa = 0.55,   
    omegaT = 0.03,  
    kappa_dlogkwHat = 0.59 ,  
    kappa_dS1m = 0.69,   
    kappa_ddS1 = 0.55,  
    rho_logkw_S1m = 0.87,  
    msigma = 0.39,  
    ssigma = 0.81)
  
  if (analyte_type == "A") {
    p=2
    beta  = c(param$beta_logkw,param$beta_S1,0,0,0,0)
    miu = c(param$logkwHat, param$S1mHat, param$dS1Hat, param$dlogkwHatA, param$dS1mHatA, param$ddS1HatA)+beta*(logP-2.2)
    msigma_sqaure = param$msigma^2
    omega = c(param$omega_logkw, param$omega_S1m, param$omega_dSa,param$kappa_dlogkwHat,param$kappa_dS1m,param$kappa_ddS1)
    rho = diag(rep(1,6))
    rho[1,2] = param$rho_logkw_S1m
    rho[2,1] = param$rho_logkw_S1m  
    tS = diag(omega)%*% rho %*%diag(omega)
  
    B = matrix(c(1,-0.4,0,0,0,0,
                 0,1,0,0,0,0,
                 0,0,1,0,0,0,
                 0,0,0,1,0,0,
                 0,0,0,0,1,0,
                 0,0,0,0,0,1),6,6,byrow = TRUE)
    
    theta= B%*%miu
    S = B%*%tS%*%t(B) 
    S_inv = inv(S)
    
    theta1 = matrix(theta[1])
    theta2 = matrix(theta[2:6])
    S11 = matrix(S[1,1])
    S12 = matrix(S[1,2:6],1,5,byrow = T)
    S21 = matrix(S[2:6,1],5,1,byrow = T)
    S22 = matrix(S[2:6,2:6],5,5,byrow = T)
    
    thetac = theta2+S21%*%inv(S11)%*%(1-theta1)
    Sc = S22 - S21%*%inv(S11)%*%S12
  } 
  
  if (analyte_type == "B") {
    p=2
    beta  = c(param$beta_logkw,param$beta_S1,0,0,0,0)
    miu = c(param$logkwHat, param$S1mHat, param$dS1Hat, param$dlogkwHatB, param$dS1mHatB, param$ddS1HatB)+beta*(logP-2.2)
    msigma_sqaure = param$msigma^2
    omega = c(param$omega_logkw, param$omega_S1m, param$omega_dSa,param$kappa_dlogkwHat,param$kappa_dS1m,param$kappa_ddS1)
    rho = diag(rep(1,6))
    rho[1,2] = param$rho_logkw_S1m
    rho[2,1] = param$rho_logkw_S1m  
    tS = diag(omega)%*% rho %*%diag(omega)
    
    B = matrix(c(1,-0.4,0,0,0,0,
                 0,1,0,0,0,0,
                 0,0,1,0,0,0,
                 0,0,0,1,0,0,
                 0,0,0,0,1,0,
                 0,0,0,0,0,1),6,6,byrow = TRUE)
    
    theta= B%*%miu
    S = B%*%tS%*%t(B) 
    S_inv = inv(S)
    
    theta1 = matrix(theta[1])
    theta2 = matrix(theta[2:6])
    S11 = matrix(S[1,1])
    S12 = matrix(S[1,2:6],1,5,byrow = T)
    S21 = matrix(S[2:6,1],5,1,byrow = T)
    S22 = matrix(S[2:6,2:6],5,5,byrow = T)
    
    thetac = theta2+S21%*%inv(S11)%*%(1-theta1)
    Sc = S22 - S21%*%inv(S11)%*%S12
  } 
  
  priors = list(
    theta = thetac,
    S = Sc)
  
  return(priors)
}

# Linear organic modifier gradient

fun_hplc_gradient_trapezoidal <- function(theta_logkw=2,
                                          theta_S1m=4,
                                          theta_dS1=0,
                                          theta_dlogkw=0,
                                          theta_dS1m=0,
                                          theta_ddS1=0,
                                          tg=20,td=1,to=1,te=0,fio=0.05,fif=0.8,
                                          m=0,xpH=0){
  td=td+0.001
  time = c(0, td, seq(td,tg+td,by=0.1))
  fi = fio-(fio-fif)/tg*(time-td)
  fi[time<td] = fio
  
  logkw = theta_logkw+xpH*theta_dlogkw											
  S1 = theta_S1m+xpH*(theta_dS1m+m*theta_ddS1)+m*theta_dS1
  S2 = 10^(-0.31+m*0.42)
  
  logki = logkw-S1*(1+S2)*fi/(1+S2*fi)
  
  ki = 10^logki 
  cumtr  = cumtrapz(time,1/to/ki)
  cumtr = as.numeric(cumtr)
  n_time = length(cumtr)
  
  if(cumtr[n_time]>=1){ # eluted before pump program ends
    trprim = suppressWarnings(interp1(cumtr,time,1))
  }else{
    trprim = (1-cumtr[n_time])* ki[n_time]*to + time[n_time]
  }
  
  tr = trprim + to + te
  return(tr)
}

fun_gra_state <- function(t, tg, td, fio, fif) {
  
  fi=fio+(fif-fio)/tg*(t-td)
  
  if (t < td) {fi = fio}
  else if (t > tg + td) {fi = fif}
  
  return(fi);
}

fun_lnki <- function(logkw, S1, S2,  fi) {
  lnki = log(10)*(logkw-S1*(1+S2)*fi/(1+S2*fi))
  return(lnki)
}  

fun_area_slope<- function(dt, lnki1, lnki2, invki1, invki2) {
  
  cki_b = c(0,0);
  
  if (invki2>1.001*invki1) {
    bo = (lnki1-lnki2)/dt
    cki = (invki2-invki1)/bo
  }
  
  else if (invki1 > 1.001 * invki2) {
    bo = (lnki2 - lnki1) / dt
    cki = (invki1 - invki2) / bo
  }
  else {
    bo = 0.001 / dt;
    cki = dt * (invki2 + invki1) / 2;
  }
  
  cki_b[1] = cki;
  cki_b[2] = bo;
  
  return(cki_b);
}


log1p_exp <- function(x) {
  indx <- .bincode(x, 
                   c(-Inf, -37, 18, 33.3, Inf), 
                   right = TRUE, 
                   include.lowest = TRUE)
  
  kk <- which(indx==1)
  if(length(kk)){x[kk]<-exp(x[kk])}
  
  kk <- which(indx==2)
  if(length(kk)){x[kk]<-log1p(exp(x[kk]))}
  
  kk <- which(indx==3)
  if(length(kk)){x[kk]<-x[kk]+exp(-x[kk])}
  
  return(x)
}

fun_hplc_gradient_steps<- function(theta_logkw=2,
                                   theta_S1m=4,
                                   theta_dS1=0,
                                   theta_dlogkw=0,
                                   theta_dS1m=0,
                                   theta_ddS1=0,
                                   tg=20,td=1,to=1,te=0,fio=0.05,fif=0.8,
                                   m=0, xpH=0) {
  
  steps = 4+6*m
  
  logkw = theta_logkw+xpH*theta_dlogkw											
  S1 = theta_S1m+xpH*(theta_dS1m+m*theta_ddS1)+m*theta_dS1
  S2 = 10^(-0.31+m*0.42)

  dt = tg / steps
  
  time1 = 0;
  time2 = td;
  
  fi1 = fun_gra_state(time1, tg, td, fio, fif);
  lnki1 = fun_lnki(logkw, S1, S2, fi1);
  lnki2=lnki1;
  
  invki1 = exp(-lnki1)/to;
  invki2 = invki1;
  
  cumki1 = 0;
  cumki2 = td*invki1;
  
  bo = 0.001/td;
  
  for (x in 1 : steps) {
    
    if (cumki2 >= 1) {
      next;
    }
    time1 = time2;
    time2 = time2 + dt;
    fi2 = fun_gra_state(time2, tg, td, fio, fif);
    lnki1 = lnki2;
    lnki2 = fun_lnki(logkw, S1, S2, fi2);
    invki1 = invki2;
    invki2 = exp(-lnki2)/to;
    cki_b = fun_area_slope(dt, lnki1, lnki2, invki1, invki2);
    cumki1 = cumki2;
    cumki2 = cumki2 + cki_b[1];
    bo = cki_b[2];
  }
  
  if (cumki2 >= 1 && cumki1==0) {
    tr = te+to+1/invki2;
  } else if (cumki2 >= 1) {
    tr = te+to+time1+log1p_exp(log((1-cumki1)*bo*to) + lnki1)/bo;
  } else if (cumki2 < 1) {
    tr = te+to+time2+(1-cumki2)/invki2;
  }
  
  return(tr);
}

compute_jacobian <- function(model, theta, X) {
  J <- jacobian(func = function(theta) model(theta, X), theta)
  return(J)
}

compute_FIM <- function(model, theta, X, S_inv) {
  J <- compute_jacobian(model, theta, X)
  FIM <- t(J) %*% J + S_inv
  return(FIM)
}

# the directional derivative 
compute_directional_derivative <- function(model, theta, X, info_matrix_inv, S_inv) {
  J <- compute_jacobian(model, theta, X)
  d <- J %*% info_matrix_inv %*% t(J)
  #d <- tr(info_matrix_inv %*% (t(J) %*% J))
  return(d)
}

fedorov_wynn_nonlinear <- function(model, theta, S_inv, p, X_candidates, chosen_indices, n_iter = 100) {
  
  n <- nrow(X_candidates) 
  X_design <- X_candidates[chosen_indices, , drop = FALSE]

  if (n_iter>0) {
  for (iter in 1:n_iter) {

    FIM = compute_FIM(model, theta, X_design, S_inv) 
    log_det_FIM <- log(det(FIM))
    
    best_gain <- 0
    best_candidate <- NULL
    best_replace_index <- NULL
    
    for (i in 1:n) {
      if (i %in% chosen_indices) next
      
      for (j in 1:p) {
        X_temp <- X_design
        X_temp[j, ] <- X_candidates[i, , drop = FALSE ]
        
        temp_info_matrix <- compute_FIM(model, theta, X_temp, S_inv) 
        temp_log_det_info_matrix <- log(det(temp_info_matrix))
      
        gain <- temp_log_det_info_matrix - log_det_FIM
        
        if (gain > best_gain) {
          best_gain <- gain
          best_candidate <- i
          best_replace_index <- j
        }
      }
    }
    
    cat("Iteration:", iter, "Best gain (ratio):", exp(best_gain), "\n")
    
    if (best_gain <= 0) {
      cat("No further improvement possible. Terminating.\n")
      break
    }
  
    chosen_indices[best_replace_index] <- best_candidate
    X_design[best_replace_index, ] <- X_candidates[best_candidate, ]
  }
  }
  FIM <- compute_FIM(model, theta, X_design, S_inv) 
  info_matrix_inv <- solve(FIM)
  det_info_matrix <- det(FIM)
  
  directional_derivatives <- NULL
  
  directional_derivatives <- sapply(1:n, function(i) {
    compute_directional_derivative(model, theta, X_candidates[i, , drop = FALSE], info_matrix_inv, S_inv)
  })

  return(list(
    chosen_indices = chosen_indices,
    X_design = X_design,
    directional_derivatives=directional_derivatives,
    det_info_matrix = det_info_matrix
  ))
}

fedorov_wynn_stan <- function(mean_log_det_stan, data, chosen_indices, n_iter = 100) {
  
  n <- nrow(data$hplcparam) 
  p=length(chosen_indices)
  
  for (iter in 1:n_iter) {
    
    log_det_FIM <- mean_log_det_stan(data,chosen_indices)
    
    best_gain <- 0
    best_candidate <- NULL
    best_replace_index <- NULL
    
    for (i in 1:n) {
      if (i %in% chosen_indices) next
      cat("n:", i, "\n")  
      for (j in 1:p) {
        chosen_indices_temp <- chosen_indices
        chosen_indices_temp[j] <- i
        
        temp_log_det_FIM <- mean_log_det_stan(data,chosen_indices_temp)
        gain <- temp_log_det_FIM - log_det_FIM
        
        if (gain > best_gain) {
          best_gain <- gain
          best_candidate <- i
          best_replace_index <- j
        }
      }
    }
    
    cat("Iteration:", iter, "Best gain (ratio):", exp(best_gain), "\n")
    
    if (best_gain <= 0) {
      cat("No further improvement possible. Terminating.\n")
      break
    }
    
    chosen_indices[best_replace_index] <- best_candidate
    cat("chosen_indices:", chosen_indices, "\n") 
  }
  

  return(list(chosen_indices = chosen_indices))
}

fedorov_wynn_stan_par <- function(mean_log_det_stan, data, init, chosen_indices, n_iter = 100) {
  
  n_rep = length(data)
  n=nrow(data[[1]]$hplcparam) 
  p=length(chosen_indices)
  if (n_iter>0) {
  for (iter in 1:n_iter) {

    log_det_FIM = parallel::mclapply(1:n_rep,\(x) mean_log_det_stan(data[[x]], init[[x]], chosen_indices, x), mc.cores=10)
    mean_log_det_FIM = mean(unlist(log_det_FIM))
    
    best_gain <- 0
    best_candidate <- NULL
    best_replace_index <- NULL
    
    for (i in 1:n) {
      if (i %in% chosen_indices) next
      cat("n:", i, "\n")  
      for (j in 1:p) {
        chosen_indices_temp <- chosen_indices
        chosen_indices_temp[j] <- i
        
        temp_log_det_FIM = parallel::mclapply(1:n_rep,\(x) mean_log_det_stan(data[[x]], init[[x]], chosen_indices_temp, x), mc.cores=10)
        temp_mean_log_det_FIM = mean(unlist(temp_log_det_FIM), na.rm = TRUE)
        
        gain <- temp_mean_log_det_FIM - mean_log_det_FIM
  
        if (gain > best_gain) {
          best_gain <- gain
          best_candidate <- i
          best_replace_index <- j
        }
      }
    }
    
    cat("Iteration:", iter, "Best average gain (ratio):", exp(best_gain), "\n")
    
    if (best_gain <= 0) {
      cat("No further improvement possible. Terminating.\n")
      break
    }
    
    chosen_indices[best_replace_index] <- best_candidate
    cat("chosen_indices:", chosen_indices, "\n") 
  }
  }
  
  log_det_FIM = parallel::mclapply(1:n_rep,\(x) mean_log_det_stan(data[[x]], init[[x]], chosen_indices, x), mc.cores=10)
  
  return(list(
    chosen_indices = chosen_indices,
    log_det_FIM = unlist(log_det_FIM)))
}

fedorov_wynn_stan_par_tgconstraint<- function(mean_log_det_stan, data, init, chosen_indices, tg_constraint=270, n_iter = 100) {
  
  
  n_rep = length(data)
  n=nrow(data[[1]]$hplcparam) 
  p=length(chosen_indices)
  if (n_iter>0) {
    for (iter in 1:n_iter) {
      
      log_det_FIM = parallel::mclapply(1:n_rep,\(x) mean_log_det_stan(data[[x]], init[[x]], chosen_indices, x), mc.cores=10)
      mean_log_det_FIM = mean(unlist(log_det_FIM))
      
      best_gain <- 0
      best_candidate <- NULL
      best_replace_index <- NULL
      
      for (i in 1:n) {
        if (i %in% chosen_indices) next
        cat("n:", i, "\n")  
        for (j in 1:p) {
          chosen_indices_temp <- chosen_indices
          chosen_indices_temp[j] <- i
          
          if (sum(data[[1]]$hplcparam[chosen_indices_temp,1])>=tg_constraint) next
          
          temp_log_det_FIM = parallel::mclapply(1:n_rep,\(x) mean_log_det_stan(data[[x]], init[[x]], chosen_indices_temp, x), mc.cores=10)
          temp_mean_log_det_FIM = mean(unlist(temp_log_det_FIM), na.rm = TRUE)
          
          gain <- temp_mean_log_det_FIM - mean_log_det_FIM
          
          if (gain > best_gain) {
            best_gain <- gain
            best_candidate <- i
            best_replace_index <- j
          }
        }
      }
      
      cat("Iteration:", iter, "Best average gain (ratio):", exp(best_gain), "\n")
      
      if (best_gain <= 0) {
        cat("No further improvement possible. Terminating.\n")
        break
      }
      
      chosen_indices[best_replace_index] <- best_candidate
      cat("chosen_indices:", chosen_indices, "\n") 
    }
  }
  
  log_det_FIM = parallel::mclapply(1:n_rep,\(x) mean_log_det_stan(data[[x]], init[[x]], chosen_indices, x), mc.cores=10)
  
  return(list(
    chosen_indices = chosen_indices,
    log_det_FIM = unlist(log_det_FIM)))
}

fun_normalized_variance= function(x){var(x)/priorsN$msigma_sqaure}
fun_cv=function(x){100*sqrt(exp(var(log(x)))-1)}
fun_lb=function(x){quantile(log10(x),0.05)}
fun_mb=function(x){quantile(log10(x),0.5)}
fun_ub=function(x){quantile(log10(x),0.95)}

fun_predict_stan_run= function(stan_model_sim, data, inits, chosen_indices, S=1) {
  
  data[[S]]$selected=chosen_indices;
  data[[S]]$nSel = length(chosen_indices);
  nAnalytes = data[[1]]$nAnalytes
  fit_opt = stan_model_sim$optimize(
    data = data[[S]],
    output_dir = here::here("stanfiles"),
    output_basename = paste0("opt-",S),
    jacobian = TRUE,
    init = function(){list(param = inits[[S]])})

  fit_results = stan_model_sim$laplace(
    data = data[[S]],
    output_dir = here::here("stanfiles"),
    output_basename = paste0("laplace-",S),
    mode = fit_opt,
    jacobian = TRUE,
    init = function(){list(param = inits[[S]])})
  
  # faster:
  draws_df= fit_results$draws(format = "df")
  
  # gra
  draws_df_subset<-draws_df[,which(colnames(draws_df) %in% grep("y_hat_sim", names(draws_df), value = TRUE))]
  draws_df_subset_var <-reshape2::melt(apply(draws_df_subset, MARGIN=2, FUN=fun_normalized_variance))
  results_gra = tibble::rownames_to_column(draws_df_subset_var, var = "rowname")%>%
    dplyr::mutate(new_id = str_extract_all(rowname, "\\d+",simplify = TRUE))%>%
    mutate(i=new_id[,1],
           c=new_id[,2])%>%
    select(i,c,variance=value)

  # izo
  draws_df_subset_M<-draws_df[,which(colnames(draws_df) %in% grep("izo_MeOH_sim", names(draws_df), value = TRUE))]
  draws_df_subset_A<-draws_df[,which(colnames(draws_df) %in% grep("izo_ACN_sim", names(draws_df), value = TRUE))]
  
  draws_df_subset_M_CV <-reshape2::melt(apply(draws_df_subset_M, MARGIN=2, FUN=fun_cv))
  draws_df_subset_A_CV <-reshape2::melt(apply(draws_df_subset_A, MARGIN=2, FUN=fun_cv))
  draws_df_subset_M_llog10tr <-reshape2::melt(apply(draws_df_subset_M, MARGIN=2, fun_lb))
  draws_df_subset_A_llog10tr <-reshape2::melt(apply(draws_df_subset_A, MARGIN=2, fun_lb))
  draws_df_subset_M_mlog10tr <-reshape2::melt(apply(draws_df_subset_M, MARGIN=2, fun_mb))
  draws_df_subset_A_mlog10tr <-reshape2::melt(apply(draws_df_subset_A, MARGIN=2, fun_mb))
  draws_df_subset_M_ulog10tr <-reshape2::melt(apply(draws_df_subset_M, MARGIN=2, FUN=fun_ub))
  draws_df_subset_A_ulog10tr <-reshape2::melt(apply(draws_df_subset_A, MARGIN=2, FUN=fun_ub))
  
  draws_df_subset_M_all = cbind(draws_df_subset_M_CV,draws_df_subset_M_llog10tr,draws_df_subset_M_ulog10tr,draws_df_subset_M_mlog10tr)
  draws_df_subset_A_all = cbind(draws_df_subset_A_CV,draws_df_subset_A_llog10tr,draws_df_subset_A_ulog10tr,draws_df_subset_A_mlog10tr)
  
  results_M = tibble::tibble(rowname = row.names(draws_df_subset_M_all), draws_df_subset_M_all, .name_repair = "universal") %>%
    dplyr::mutate(new_id = str_extract_all(rowname, "\\d+",simplify = TRUE)) %>%
    mutate(i=new_id[,1],
           c=new_id[,2])%>%
    select(i,c,CV=`value...2`,llog10tr=`value...3`,  ulog10tr=`value...4`, mlog10tr=`value...5`) %>%
    mutate(mod = "MeOH")
  
  results_A = tibble::tibble(rowname = row.names(draws_df_subset_A_all), draws_df_subset_A_all, .name_repair = "universal") %>%
    dplyr::mutate(new_id = str_extract_all(rowname, "\\d+",simplify = TRUE)) %>%
    mutate(i=new_id[,1],
           c=new_id[,2])%>%
    select(i,c,CV=`value...2`,llog10tr=`value...3`,  ulog10tr=`value...4`, mlog10tr=`value...5`) %>%
    mutate(mod = "ACN")

  results_izo = rbind(results_M,results_A)
  return(list(results_gra = results_gra,
              results_izo = results_izo))
}

fun_predict_stan_run_pH= function(stan_model_sim_pH, data_pH, inits_pH, chosen_indices, S=1) {
  
  data_pH[[S]]$selected=chosen_indices;
  data_pH[[S]]$nSel = length(chosen_indices);
  nAnalytes = data_pH[[1]]$nAnalytes
  fit_opt = stan_model_sim_pH$optimize(
    data = data_pH[[S]],
    output_dir = here::here("stanfiles"),
    output_basename = paste0("opt-",S),
    jacobian = TRUE,
    init = function(){list(param = inits_pH[[S]])})
  
  fit_results = stan_model_sim_pH$laplace(
    data = data_pH[[S]],
    output_dir = here::here("stanfiles"),,
    output_basename = paste0("laplace-",S),
    mode = fit_opt,
    jacobian = TRUE,
    init = function(){list(param = inits_pH[[S]])})
  
  Sys.sleep(10)
  
  draws_df= fit_results$draws(format = "df")
  
  # gra
  draws_df_subset<-draws_df[,which(colnames(draws_df) %in% grep("y_hat_sim", names(draws_df), value = TRUE))]
  draws_df_subset_var <-reshape2::melt(apply(draws_df_subset, MARGIN=2, FUN=fun_normalized_variance))
  results_gra = tibble::rownames_to_column(draws_df_subset_var, var = "rowname")%>%
    dplyr::mutate(new_id = stringr::str_extract_all(rowname, "\\d+",simplify = TRUE))%>%
    mutate(i=new_id[,1],
           c=new_id[,2])%>%
    select(i,c,variance=value)
  
  # izo
  draws_df_subset_MN<-draws_df[,which(colnames(draws_df) %in% grep("izo_MeOH_sim_N", names(draws_df), value = TRUE))]
  draws_df_subset_AN<-draws_df[,which(colnames(draws_df) %in% grep("izo_ACN_sim_N", names(draws_df), value = TRUE))]
  draws_df_subset_MI<-draws_df[,which(colnames(draws_df) %in% grep("izo_MeOH_sim_I", names(draws_df), value = TRUE))]
  draws_df_subset_AI<-draws_df[,which(colnames(draws_df) %in% grep("izo_ACN_sim_I", names(draws_df), value = TRUE))]
  
  
  draws_df_subset_MN_CV <-reshape2::melt(apply(draws_df_subset_MN, MARGIN=2, FUN=fun_cv))
  draws_df_subset_AN_CV <-reshape2::melt(apply(draws_df_subset_AN, MARGIN=2, FUN=fun_cv))
  draws_df_subset_MN_llog10tr <-reshape2::melt(apply(draws_df_subset_MN, MARGIN=2, FUN=fun_lb))
  draws_df_subset_AN_llog10tr <-reshape2::melt(apply(draws_df_subset_AN, MARGIN=2, FUN=fun_lb))
  draws_df_subset_MN_mlog10tr <-reshape2::melt(apply(draws_df_subset_MN, MARGIN=2, FUN=fun_mb))
  draws_df_subset_AN_mlog10tr <-reshape2::melt(apply(draws_df_subset_AN, MARGIN=2, FUN=fun_mb))
  draws_df_subset_MN_ulog10tr <-reshape2::melt(apply(draws_df_subset_MN, MARGIN=2, FUN=fun_ub))
  draws_df_subset_AN_ulog10tr <-reshape2::melt(apply(draws_df_subset_AN, MARGIN=2, FUN=fun_ub))
  
  draws_df_subset_MI_CV <-reshape2::melt(apply(draws_df_subset_MI, MARGIN=2, FUN=fun_cv))
  draws_df_subset_AI_CV <-reshape2::melt(apply(draws_df_subset_AI, MARGIN=2, FUN=fun_cv))
  draws_df_subset_MI_llog10tr <-reshape2::melt(apply(draws_df_subset_MI, MARGIN=2, FUN=fun_lb))
  draws_df_subset_AI_llog10tr <-reshape2::melt(apply(draws_df_subset_AI, MARGIN=2, FUN=fun_lb))
  draws_df_subset_MI_mlog10tr <-reshape2::melt(apply(draws_df_subset_MI, MARGIN=2, FUN=fun_mb))
  draws_df_subset_AI_mlog10tr <-reshape2::melt(apply(draws_df_subset_AI, MARGIN=2, FUN=fun_mb))
  draws_df_subset_MI_ulog10tr <-reshape2::melt(apply(draws_df_subset_MI, MARGIN=2, FUN=fun_ub))
  draws_df_subset_AI_ulog10tr <-reshape2::melt(apply(draws_df_subset_AI, MARGIN=2, FUN=fun_ub))
  
  draws_df_subset_MN_all = cbind(draws_df_subset_MN_CV,draws_df_subset_MN_llog10tr,draws_df_subset_MN_ulog10tr,draws_df_subset_MN_mlog10tr)
  draws_df_subset_MI_all = cbind(draws_df_subset_MI_CV,draws_df_subset_MI_llog10tr,draws_df_subset_MI_ulog10tr,draws_df_subset_MI_mlog10tr)
  draws_df_subset_AN_all = cbind(draws_df_subset_AN_CV,draws_df_subset_AN_llog10tr,draws_df_subset_AN_ulog10tr,draws_df_subset_AN_mlog10tr)
  draws_df_subset_AI_all = cbind(draws_df_subset_AI_CV,draws_df_subset_AI_llog10tr,draws_df_subset_AI_ulog10tr,draws_df_subset_AI_mlog10tr)
  
  results_MN= tibble::tibble(rowname = row.names(draws_df_subset_MN_all), draws_df_subset_MN_all, .name_repair = "universal") %>%
    dplyr::mutate(new_id = str_extract_all(rowname, "\\d+",simplify = TRUE)) %>%
    mutate(i=new_id[,1],
           c=new_id[,2])%>%
    select(i,c,CV=`value...2`,llog10tr=`value...3`,  ulog10tr=`value...4`, mlog10tr=`value...5`) %>%
    mutate(mod = "MeOH")%>%
    mutate(form = "N")
  
  results_MI= tibble::tibble(rowname = row.names(draws_df_subset_MI_all), draws_df_subset_MI_all, .name_repair = "universal") %>%
    dplyr::mutate(new_id = str_extract_all(rowname, "\\d+",simplify = TRUE)) %>%
    mutate(i=new_id[,1],
           c=new_id[,2])%>%
    select(i,c,CV=`value...2`,llog10tr=`value...3`, ulog10tr=`value...4`, mlog10tr=`value...5`) %>%
    mutate(mod = "MeOH")%>%
    mutate(form = "I")
  
  results_AN = tibble::tibble(rowname = row.names(draws_df_subset_AN_all), draws_df_subset_AN_all, .name_repair = "universal") %>%
    dplyr::mutate(new_id = str_extract_all(rowname, "\\d+",simplify = TRUE)) %>%
    mutate(i=new_id[,1],
           c=new_id[,2])%>%
    select(i,c,CV=`value...2`,llog10tr=`value...3`,  ulog10tr=`value...4`, mlog10tr=`value...5`) %>%
    mutate(mod = "ACN")%>%
    mutate(form = "N")
  
  results_AI = tibble::tibble(rowname = row.names(draws_df_subset_AI_all), draws_df_subset_AI_all, .name_repair = "universal") %>%
    dplyr::mutate(new_id = str_extract_all(rowname, "\\d+",simplify = TRUE)) %>%
    mutate(i=new_id[,1],
           c=new_id[,2])%>%
    select(i,c,CV=`value...2`,llog10tr=`value...3`,  ulog10tr=`value...4`, mlog10tr=`value...5`) %>%
    mutate(mod = "ACN")%>%
    mutate(form = "I")
  
  results_izo = rbind(results_MN,results_AN,results_MI,results_AI)
  
  return(list(results_gra = results_gra,
              results_izo = results_izo))
}

fun_predict_stan_run_MeOH= function(stan_model_sim_MeOH, data_MeOH, inits_MeOH, chosen_indices, S=1) {
  
  data_MeOH[[S]]$selected=chosen_indices;
  data_MeOH[[S]]$nSel = length(chosen_indices);
  nAnalytes = data_MeOH[[1]]$nAnalytes
  fit_opt = stan_model_sim_MeOH$optimize(
    data = data_MeOH[[S]],
    output_dir = here::here("stanfiles"),
    output_basename = paste0("opt-",S),
    jacobian = TRUE,
    init = function(){list(param = inits_MeOH[[S]])})
  
  fit_results = stan_model_sim_MeOH$laplace(
    data = data_MeOH[[S]],
    output_dir = here::here("stanfiles"),
    output_basename = paste0("laplace-",S),
    mode = fit_opt,
    jacobian = TRUE,
    init = function(){list(param = inits_MeOH[[S]])})

  Sys.sleep(10)
  
  draws_df= fit_results$draws(format = "df")
  
  # gra
  draws_df_subset<-draws_df[,which(colnames(draws_df) %in% grep("y_hat_sim", names(draws_df), value = TRUE))]
  draws_df_subset_var <-reshape2::melt(apply(draws_df_subset, MARGIN=2, FUN=fun_normalized_variance))
  results_gra = tibble::rownames_to_column(draws_df_subset_var, var = "rowname")%>%
    dplyr::mutate(new_id = str_extract_all(rowname, "\\d+",simplify = TRUE))%>%
    mutate(i=new_id[,1],
           c=new_id[,2])%>%
    select(i,c,variance=value)
  
  # izo
  draws_df_subset_M<-draws_df[,which(colnames(draws_df) %in% grep("izo_MeOH_sim", names(draws_df), value = TRUE))]

  draws_df_subset_M_CV <-reshape2::melt(apply(draws_df_subset_M, MARGIN=2, FUN=fun_cv))
  draws_df_subset_M_llog10tr <-reshape2::melt(apply(draws_df_subset_M, MARGIN=2, FUN=fun_lb))
  draws_df_subset_M_mlog10tr <-reshape2::melt(apply(draws_df_subset_M, MARGIN=2, FUN=fun_mb))
  draws_df_subset_M_ulog10tr <-reshape2::melt(apply(draws_df_subset_M, MARGIN=2, FUN=fun_ub))
  
  draws_df_subset_M_all = cbind(draws_df_subset_M_CV,draws_df_subset_M_llog10tr,draws_df_subset_M_ulog10tr,draws_df_subset_M_mlog10tr)

  results_M = tibble::tibble(rowname = row.names(draws_df_subset_M_all), draws_df_subset_M_all, .name_repair = "universal") %>%
    dplyr::mutate(new_id = str_extract_all(rowname, "\\d+",simplify = TRUE)) %>%
    mutate(i=new_id[,1],
           c=new_id[,2])%>%
    select(i,c,CV=`value...2`,llog10tr=`value...3`,  ulog10tr=`value...4`, mlog10tr=`value...5`) %>%
    mutate(mod = "MeOH")

  
  results_izo = rbind(results_M)
  return(list(results_gra = results_gra,
              results_izo = results_izo))
}

fun_predict_stan_run_ACN= function(stan_model_sim_ACN, data_ACN, inits_ACN, chosen_indices, S=1) {
  
  data_ACN[[S]]$selected=chosen_indices;
  data_ACN[[S]]$nSel = length(chosen_indices);
  nAnalytes = data_ACN[[1]]$nAnalytes
  fit_opt = stan_model_sim_ACN$optimize(
    data = data_ACN[[S]],
    output_dir = here::here("stanfiles"),
    output_basename = paste0("opt-",S),
    jacobian = TRUE,
    init = function(){list(param = inits_ACN[[S]])})
  
  fit_results = stan_model_sim_ACN$laplace(
    data = data_ACN[[S]],
    output_dir = here::here("stanfiles"),
    output_basename = paste0("laplace-",S),
    mode = fit_opt,
    jacobian = TRUE,
    init = function(){list(param = inits_ACN[[S]])})
  
  Sys.sleep(10)
  
  draws_df= fit_results$draws(format = "df")
  
  # gra
  draws_df_subset<-draws_df[,which(colnames(draws_df) %in% grep("y_hat_sim", names(draws_df), value = TRUE))]
  draws_df_subset_var <-reshape2::melt(apply(draws_df_subset, MARGIN=2, FUN=fun_normalized_variance))
  results_gra = tibble::rownames_to_column(draws_df_subset_var, var = "rowname")%>%
    dplyr::mutate(new_id = str_extract_all(rowname, "\\d+",simplify = TRUE))%>%
    mutate(i=new_id[,1],
           c=new_id[,2])%>%
    select(i,c,variance=value)
  
  # izo
  draws_df_subset_A<-draws_df[,which(colnames(draws_df) %in% grep("izo_ACN_sim", names(draws_df), value = TRUE))]
  
  draws_df_subset_A_CV <-reshape2::melt(apply(draws_df_subset_A, MARGIN=2, FUN=fun_cv))
  draws_df_subset_A_llog10tr <-reshape2::melt(apply(draws_df_subset_A, MARGIN=2, FUN=fun_lb))
  draws_df_subset_A_mlog10tr <-reshape2::melt(apply(draws_df_subset_A, MARGIN=2, FUN=fun_mb))
  draws_df_subset_A_ulog10tr <-reshape2::melt(apply(draws_df_subset_A, MARGIN=2, FUN=fun_ub))
  
  draws_df_subset_A_all = cbind(draws_df_subset_A_CV,draws_df_subset_A_llog10tr,draws_df_subset_A_ulog10tr,draws_df_subset_A_mlog10tr)

  results_A = tibble::tibble(rowname = row.names(draws_df_subset_A_all), draws_df_subset_A_all, .name_repair = "universal") %>%
  dplyr::mutate(new_id = str_extract_all(rowname, "\\d+",simplify = TRUE)) %>%
  mutate(i=new_id[,1],
         c=new_id[,2])%>%
  select(i,c,CV=`value...2`,llog10tr=`value...3`,  ulog10tr=`value...4`, mlog10tr=`value...5`) %>%
  mutate(mod = "ACN")
  
  results_izo = rbind(results_A)
  return(list(results_gra = results_gra,
              results_izo = results_izo))
}

fun_predict_stan_results = function(stan_model_sim, data, inits, chosen_indices, n_rep) {
  
  nAnalytes = data[[1]]$nAnalytes
  
  results_temp= bettermc::mclapply(1:n_rep,\(i) fun_predict_stan_run(stan_model_sim, 
                                                                     data, 
                                                                     inits,
                                                                     chosen_indices, 
                                                                     S=i), mc.cores=10)
  
  results_gra <- lapply(results_temp, `[[`, "results_gra")
  
  results=bind_rows(results_gra, .id = "S")%>%
    unnest(cols = c()) %>%
    mutate(i=as.factor(as.numeric(i)+(as.numeric(S)-1)*nAnalytes))
  
  iresults=results %>%
    group_by(i) %>%
    summarise(ivar = max(variance),
              imvar = mean(variance))
  
  cresults = results %>%
    group_by(c) %>%
    summarise(qvar = median(variance))
  
  results_izo <- lapply(results_temp, `[[`, "results_izo")
  
  results=bind_rows(results_izo, .id = "S")%>%
    unnest(cols = c()) %>%
    mutate(i=as.factor(as.numeric(i)+(as.numeric(S)-1)*nAnalytes))
  
  return(list(iresults=iresults,
              cresults=cresults,
              izoresults=results))
} 

fun_predict_stan_results_pH = function(stan_model_sim_pH, data_pH, inits_pH, chosen_indices, n_rep) {
  
  nAnalytes = data_pH[[1]]$nAnalytes
  
  results_temp= bettermc::mclapply(1:n_rep,\(i) fun_predict_stan_run_pH(stan_model_sim_pH, 
                                                                        data_pH, 
                                                                     inits_pH,
                                                                     chosen_indices, 
                                                                     S=i), mc.cores=10)
  

  results_gra <- lapply(results_temp, `[[`, "results_gra")
  
  results=bind_rows(results_gra, .id = "S")%>%
    unnest(cols = c()) %>%
    mutate(i=as.factor(as.numeric(i)+(as.numeric(S)-1)*nAnalytes))
  
  iresults=results %>%
    group_by(i) %>%
    summarise(ivar = max(variance))
  
  cresults = results %>%
    group_by(c) %>%
    summarise(qvar = median(variance))
  
  results_izo <- lapply(results_temp, `[[`, "results_izo")
  
  results=bind_rows(results_izo, .id = "S")%>%
    unnest(cols = c()) %>%
    mutate(i=as.factor(as.numeric(i)+(as.numeric(S)-1)*nAnalytes))
  
  return(list(iresults=iresults,
              cresults=cresults,
              izoresults=results))
} 

fun_predict_stan_results_MeOH = function(stan_model_sim_MeOH, data_MeOH, inits_MeOH, chosen_indices, n_rep) {
  
  nAnalytes = data_MeOH[[1]]$nAnalytes
  
  results_temp= bettermc::mclapply(1:n_rep,\(i) fun_predict_stan_run_MeOH(stan_model_sim_MeOH, 
                                                                        data_MeOH, 
                                                                        inits_MeOH,
                                                                        chosen_indices, 
                                                                        S=i), mc.cores=10)
  
  results_gra <- lapply(results_temp, `[[`, "results_gra")
  
  results=bind_rows(results_gra, .id = "S")%>%
    unnest(cols = c()) %>%
    mutate(i=as.factor(as.numeric(i)+(as.numeric(S)-1)*nAnalytes))
  
  iresults=results %>%
    group_by(i) %>%
    summarise(ivar = max(variance))
  
  cresults = results %>%
    group_by(c) %>%
    summarise(qvar = median(variance))
  
  results_izo <- lapply(results_temp, `[[`, "results_izo")
  
  results=bind_rows(results_izo, .id = "S")%>%
    unnest(cols = c()) %>%
    mutate(i=as.factor(as.numeric(i)+(as.numeric(S)-1)*nAnalytes))
  
  return(list(iresults=iresults,
              cresults=cresults,
              izoresults=results))
}           

fun_predict_stan_results_ACN = function(stan_model_sim_ACN, data_ACN, inits_ACN, chosen_indices, n_rep) {
  
  nAnalytes = data_ACN[[1]]$nAnalytes
  
  results_temp= bettermc::mclapply(1:n_rep,\(i) fun_predict_stan_run_ACN(stan_model_sim_ACN, 
                                                                          data_ACN, 
                                                                          inits_ACN,
                                                                          chosen_indices, 
                                                                          S=i), mc.cores=10)
  
  results_gra <- lapply(results_temp, `[[`, "results_gra")
  
  results=bind_rows(results_gra, .id = "S")%>%
    unnest(cols = c()) %>%
    mutate(i=as.factor(as.numeric(i)+(as.numeric(S)-1)*nAnalytes))
  
  iresults=results %>%
    group_by(i) %>%
    summarise(ivar = max(variance))
  
  cresults = results %>%
    group_by(c) %>%
    summarise(qvar = median(variance))
  
  results_izo <- lapply(results_temp, `[[`, "results_izo")
  
  results=bind_rows(results_izo, .id = "S")%>%
    unnest(cols = c()) %>%
    mutate(i=as.factor(as.numeric(i)+(as.numeric(S)-1)*nAnalytes))
  
  return(list(iresults=iresults,
              cresults=cresults,
              izoresults=results))
}  