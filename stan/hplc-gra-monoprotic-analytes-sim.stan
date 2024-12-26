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
  int nAnalytes;  // number of analytes
  int nObs;// number of observations
  array[nObs] int<lower=1> steps;
  array[nObs] vector[8] hplcparam; // [tg, td, to, te, fio, fik, mod, pH]
  int<lower=0, upper=1> run_estimation; // 0 for prior predictive, 1 for estimation
  real<lower=0, upper=1> include_prior; // small number of no prior.
  int nSel;
  array[nSel] int<lower=1> selected;
  matrix[nAnalytes,nObs] trobs; 
  vector[nAnalytes] logPobs; 
  array[nAnalytes] int A1B2;
}

transformed data{
matrix[nAnalytes,nSel] trobs_sel = trobs[1:nAnalytes, selected];
vector[11] fi_sim =[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]';
array[2] cov_matrix[6] Omega;
matrix[2,6] prior_theta;
prior_theta = [[3.60,4.92,0.61,-0.79,0.17,0.28],[3.60,4.92,0.61,-0.97,0.12,-0.67]];
Omega[1] = 1/include_prior*[[0.846400,0.744372,0.0000,0.0000,0.0000,0.0000],
[0.744372,0.864900,0.0000,0.0000,0.0000,0.0000],
[0.000000,0.000000,0.3025,0.0000,0.0000,0.0000],
[0.000000,0.000000,0.0000,0.3481,0.0000,0.0000],
[0.000000,0.000000,0.0000,0.0000,0.4761,0.0000],
[0.000000,0.000000,0.0000,0.0000,0.0000,0.3025]];
Omega[2] = 1/include_prior*[[0.846400,0.744372,0.0000,0.0000,0.0000,0.0000],
[0.744372,0.864900,0.0000,0.0000,0.0000,0.0000],
[0.000000,0.000000,0.3025,0.0000,0.0000,0.0000],
[0.000000,0.000000,0.0000,0.3481,0.0000,0.0000],
[0.000000,0.000000,0.0000,0.0000,0.4761,0.0000],
[0.000000,0.000000,0.0000,0.0000,0.0000,0.3025]];
real<lower=0> sigma = 0.39;
}

parameters {
array[nAnalytes] vector[6] param; // logkw, S1m, dS1
}

transformed parameters {
array[nAnalytes] vector[6] miu;
matrix[nAnalytes, nSel] y_hat;
   
   
for (i in 1 : nAnalytes) {
  miu[i,1] = prior_theta[A1B2[i],1] + 0.84*(logPobs[i]-2.2);
  miu[i,2] = prior_theta[A1B2[i],2] + 0.51*(logPobs[i]-2.2);
  miu[i,3] = prior_theta[A1B2[i],3];
  miu[i,4] = prior_theta[A1B2[i],4];
  miu[i,5] = prior_theta[A1B2[i],5];
  miu[i,6] = prior_theta[A1B2[i],6];
 }
 
 {
 for (z in 1 : nSel){
   for (i in 1 : nAnalytes){
    real logkw = param[i,1]	+ hplcparam[selected[z],8]*param[i,4];								
    real S1 = param[i,2] + hplcparam[selected[z],8]*(param[i,5]+hplcparam[selected[z],7]*param[i,6])+hplcparam[selected[z],7]*param[i,3];
    real S2 = 10^(-0.31+hplcparam[selected[z],7]*0.42);
    y_hat[i,z] = chromgratrapz(steps[selected[z]], logkw, S1, S2, hplcparam[selected[z]]);
  }
 }
  }
}


model {

 for (i in 1 : nAnalytes) {
  param[i] ~ multi_normal(miu[i], Omega[A1B2[i]]); 
  }
  
 if (run_estimation == 1) {
 to_vector(trobs_sel) ~ normal(to_vector(y_hat), sigma);
 }
}

generated quantities {
   
array[nAnalytes] vector[nObs] y_hat_sim;
array[nAnalytes] vector[11] izo_MeOH_sim_N;
array[nAnalytes] vector[11] izo_MeOH_sim_I;
array[nAnalytes] vector[11] izo_ACN_sim_N;
array[nAnalytes] vector[11] izo_ACN_sim_I;
 
 {
 for (z in 1 : nObs){
   for (i in 1 : nAnalytes){
    real logkw_sim = param[i,1]	+ hplcparam[z,8]*param[i,4];								
    real S1_sim = param[i,2] + hplcparam[z,8]*(param[i,5]+hplcparam[z,7]*param[i,6])+hplcparam[z,7]*param[i,3];
    real S2_sim = 10^(-0.31+hplcparam[z,7]*0.42);
    y_hat_sim[i,z] = chromgratrapz(steps[z], logkw_sim, S1_sim, S2_sim, hplcparam[z]);

  }
 }

for (i in 1 : nAnalytes){   
 for (z in 1 : 11){
   real logkw_sim = param[i,1];
   real S1_sim = param[i,2] ;
   real S2_sim = 10^(-0.31);
   izo_MeOH_sim_N[i,z] = hplcparam[z,4] + hplcparam[z,3] *(1+exp(funlnki(logkw_sim, S1_sim, S2_sim, fi_sim[z])));
  }
}
  
for (i in 1 : nAnalytes){   
 for (z in 1 : 11){
   real logkw_sim = param[i,1];
   real S1_sim = param[i,2] + param[i,3];
   real S2_sim = 10^(-0.31 + 0.42);
   izo_ACN_sim_N[i,z] = hplcparam[z,4] + hplcparam[z,3] *(1+exp(funlnki(logkw_sim, S1_sim, S2_sim, fi_sim[z])));
  }
}

for (i in 1 : nAnalytes){   
 for (z in 1 : 11){
   real logkw_sim = param[i,1] + param[i,4];
   real S1_sim = param[i,2] + param[i,5];
   real S2_sim = 10^(-0.31);
   izo_MeOH_sim_I[i,z] = hplcparam[z,4] + hplcparam[z,3] *(1+exp(funlnki(logkw_sim, S1_sim, S2_sim, fi_sim[z])));
  }
}

for (i in 1 : nAnalytes){   
 for (z in 1 : 11){
   real logkw_sim = param[i,1] + param[i,4];
   real  S1_sim = param[i,2] + param[i,3] + param[i,5] + param[i,6];
   real S2_sim = 10^(-0.31 + 0.42);
   izo_ACN_sim_I[i,z] = hplcparam[z,4] + hplcparam[z,3] *(1+exp(funlnki(logkw_sim, S1_sim, S2_sim, fi_sim[z])));
  }
}
}
}
