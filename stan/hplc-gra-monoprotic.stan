functions {
  
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
  int nObs;// number of observations
  array[nObs] int<lower=1> steps;
  array[nObs] vector[8] hplcparam; // [tg, td, to, te, fio, fik, mod, xpH]
  real logPobs;
  cov_matrix[6] Omega;
  real<lower=0> sigma;
  vector[nObs] trobs; 
  int<lower=0, upper=1> run_estimation; // 0 for prior predictive, 1 for estimation
  
  // simulations
  int nObs_sim;// number of observations
  array[nObs_sim] int<lower=1> steps_sim;
  array[nObs_sim] vector[8] hplcparam_sim; // [tg, td, to, te, fio, fik, mod, xpH]
}

transformed data {
vector[11] fi_sim =[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]';
}

parameters {
  vector[6] param; // logkw, S1m, dS1
}

transformed parameters {

  vector[6] miu;
  vector[nObs] y_hat;
  
  miu[1] = 3.60 + 0.84 * (logPobs-2.2);
  miu[2] = 4.92 + 0.51 * (logPobs-2.2);
  miu[3] = 0.61;
  miu[4] =-0.79;
  miu[5] = 0.17;
  miu[6] = 0.28;

 {
   real logkw;
   real S1;
   real S2; 
   
 for (z in 1 : nObs){
     logkw = param[1]	+ hplcparam[z,8]*param[4];									
     S1 = param[2] + hplcparam[z,8]*(param[5]+hplcparam[z,7]*param[6])+hplcparam[z,7]*param[3];
     S2 = 10^(-0.31+hplcparam[z,7]*0.42);
  
   y_hat[z] = chromgratrapz(steps[z], logkw, S1, S2, hplcparam[z]);
  }
  }
  
}


model {
  
param ~ multi_normal(miu, Omega); 

 if (run_estimation == 1) {
 trobs ~ normal(y_hat, sigma);
 }
}

generated quantities {
   
 vector[nObs_sim] y_hat_sim;
 vector[11] izo_MeOH_sim_N;
 vector[11] izo_MeOH_sim_I;
 vector[11] izo_ACN_sim_N;
 vector[11] izo_ACN_sim_I;

 {
 for (z in 1 : nObs_sim){
     real logkw_sim = param[1]	+ hplcparam_sim[z,8]*param[4];								
     real S1_sim = param[2] + hplcparam_sim[z,8]*(param[5]+hplcparam_sim[z,7]*param[6])+hplcparam_sim[z,7]*param[3];
     real S2_sim = 10^(-0.31+hplcparam_sim[z,7]*0.42);
   y_hat_sim[z] = chromgratrapz(steps_sim[z], logkw_sim, S1_sim, S2_sim, hplcparam_sim[z]);
  }

for (z in 1 : 11){
   real logkw_sim = param[1];
   real S1_sim = param[2] ;
   real S2_sim = 10^(-0.31);
   izo_MeOH_sim_N[z] = hplcparam_sim[z,4] + hplcparam_sim[z,3] *(1+exp(funlnki(logkw_sim, S1_sim, S2_sim, fi_sim[z])));
  }
 for (z in 1 : 11){
   real logkw_sim = param[1];
   real S1_sim = param[2] + param[3];
   real S2_sim = 10^(-0.31 + 0.42);
   izo_ACN_sim_N[z] = hplcparam_sim[z,4] + hplcparam_sim[z,3] *(1+exp(funlnki(logkw_sim, S1_sim, S2_sim, fi_sim[z])));
  }
  for (z in 1 : 11){
   real logkw_sim = param[1] + param[4];
   real    S1_sim = param[2] + param[5];
   real S2_sim = 10^(-0.31);
   izo_MeOH_sim_I[z] = hplcparam_sim[z,4] + hplcparam_sim[z,3] *(1+exp(funlnki(logkw_sim, S1_sim, S2_sim, fi_sim[z])));
  }
 for (z in 1 : 11){
   real logkw_sim = param[1] + param[4];
   real  S1_sim = param[2] + param[3] + param[5] + param[6];
   real S2_sim = 10^(-0.31 + 0.42);
   izo_ACN_sim_I[z] = hplcparam_sim[z,4] + hplcparam_sim[z,3] *(1+exp(funlnki(logkw_sim, S1_sim, S2_sim, fi_sim[z])));
  } 
}
}
