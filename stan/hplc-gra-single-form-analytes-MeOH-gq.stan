

functions {
  
real normal_lub_rng(real mu, real sigma, real lb, real ub) {
  real p_lb = normal_cdf(lb | mu, sigma);
  real p_ub = normal_cdf(ub | mu, sigma);
  real u = uniform_rng(p_lb, p_ub);
  real y = mu + sigma * inv_Phi(u);
  return y;
}

  // fi at a given time at column inlet
  real gra_state(real t, vector hplcparam) {
    real tg = hplcparam[1];
    real td = hplcparam[2];
    real fio = hplcparam[5];
    real fik = hplcparam[6];
    real fi;
    
    fi=fio+(fik-fio)/tg*(t-td);
    
    if (t < td) {
      fi = fio;
    } else if (t > tg + td) {
      fi = fik;
    }
    return fi;
  }
  
 real funlnki(real logkw, real S1, real S2, real fi) {
    real lnki;
    lnki = log(10) *(logkw-S1*(1+S2)*fi/(1+S2*fi));
    return lnki;
  }  
 
  vector areaandslope(real dt, real lnki1, real lnki2, real invki1, real invki2) {
    vector[2] cki_b;
    real bo;
    real cki;
    
    if (invki2 > 1.001 * invki1) {
      bo = (lnki1 - lnki2) / dt;
      cki = (invki2 - invki1) / bo;
    }
    
    else if (invki1 > 1.001 * invki2) {
      bo = (lnki2 - lnki1) / dt;
      cki = (invki1 - invki2) / bo;
    }
    else {
      bo = 0.001 / dt;
      cki = dt * (invki2 + invki1) / 2;
    }
    
    cki_b[1] = cki;
    cki_b[2] = bo;
    
    return cki_b;
  }
  
  real chromgratrapz(int steps, real logkw, real S1,  real S2, vector hplcparam) {
                       
    real tg = hplcparam[1];
    real td = hplcparam[2];
    real to = hplcparam[3];
    real te = hplcparam[4];
    
    real time1;
    real time2;
    real fi1;
    real fi2;
    real lnki1;
    real lnki2;
    real invki1;
    real invki2;
    vector[2] cki_b;
    real cumki1;
    real cumki2;
    real bo;
    real tr;
    real dt;
    
    dt = tg / steps;
    
    time1 = 0;
    time2 = td;
    
    fi1 = gra_state(time1, hplcparam);
    lnki1 = funlnki(logkw, S1, S2, fi1);
    lnki2=lnki1;
    
    invki1 = exp(-lnki1)/to;
    invki2 = invki1;
    
    cumki1 = 0;
    cumki2 = td*invki1; 
    
    bo = 0.001 / td;
    
    for (x in 1 : steps) {
      if (cumki2 >= 1) {
        continue;
      }
      time1 = time2;
      time2 += dt;
      fi2 = gra_state(time2, hplcparam);
      lnki1 = lnki2;
      lnki2 = funlnki(logkw, S1, S2, fi2);
      invki1 = invki2;
      invki2 = exp(-lnki2)/to;
      cki_b = areaandslope(dt, lnki1, lnki2, invki1, invki2);
      cumki1 = cumki2;
      cumki2 += cki_b[1];
      bo = cki_b[2];
    }
    
    if (cumki2 >= 1 && cumki1==0) {
      tr = te+to+1/invki2;
    } else if (cumki2 >= 1) {
      tr = te+to+time1+log1p_exp(log((1-cumki1)*bo*to) + lnki1)/bo;
    } else if (cumki2 < 1) {
      tr = te+to+time2+(1-cumki2)/invki2;
    }
    return tr;
  }
}

data {
  int nAnalytes;  // number of analytes
  int nObs;// number of observations
  array[nObs] int<lower=1> steps;
  array[nObs] vector[7] hplcparam; // [tg, td, to, te, fio, fik, mod]
  real<lower=0, upper=1> include_prior; // small number of no prior.

  int<lower=0, upper=1> run_estimation; // 0 for prior predictive, 1 for estimation
  int nSel;
  array[nSel] int<lower=1> selected;
}

transformed data {
array[3] cov_matrix[2] Omega;
vector[3] miu_log = [2.2,2.2,2.2]';
matrix[3,2] prior_theta;
prior_theta = [[3.60,4.92],[2.81,5.09],[2.63,5.04]];
Omega[1] = 1/include_prior*[[0.846400, 0.744372], [0.744372, 0.864900]];
Omega[2] = 1/include_prior*[[1.194500, 0.744372], [0.744372, 1.341000]];
Omega[3] = 1/include_prior*[[1.194500, 0.744372], [0.744372, 1.341000]];
real<lower=0> sigma = 0.39;
}

parameters {
}

transformed parameters {
}

model {
}

generated quantities {
array[nAnalytes] vector[2] miu;
matrix[nAnalytes,nSel] y_hat_sim;
matrix[nAnalytes,nSel] trobs; 
vector[nAnalytes] logPobs; 
array[nAnalytes] int N1A2B3;
array[nAnalytes] vector[2] paramsim; // logkw, S1m, dS1
 
 for (i in 1 : nAnalytes) {
  N1A2B3[i] = categorical_rng([0.4,0.3,0.3]');
  logPobs[i] = normal_lub_rng(miu_log[N1A2B3[i]],2,0.5,7);
  miu[i,1] = prior_theta[N1A2B3[i],1] + 0.84*(logPobs[i]-2.2);
  miu[i,2] = prior_theta[N1A2B3[i],2] + 0.51*(logPobs[i]-2.2);
  paramsim[i] = multi_normal_rng(miu[i], Omega[N1A2B3[i]]);
 }
 
 {
   real logkws;
   real S1s;
   real S2s; 
 for (z in 1 : nSel){
   for (i in 1 : nAnalytes){
   logkws = paramsim[i,1];
   S1s = paramsim[i,2];
   S2s = 10^(-0.31);
   y_hat_sim[i,z] = chromgratrapz(steps[selected[z]], logkws, S1s, S2s, hplcparam[selected[z]]);
   trobs[i,z] = normal_rng( y_hat_sim[i,z], sigma);
  }
 }
}
}
