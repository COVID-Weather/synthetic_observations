data {
  real<lower=0> POP;          	
  int<lower=0>  N_obs;         
  int<lower=0>  y[N_obs]; 
  real          E0;
  real<lower=0> Ntest[N_obs];
  //real s;
}
parameters {
  real<lower=0.0> beta[N_obs];
  //real<lower=0.0> beta;
  //real<lower=0.0> gamma;
  //real<lower=0.0> sigma;
  real<lower=0.0,upper=1.0> s;
  //real<lower=0.0> E0;
  //real<lower=1E-15,upper=0.1> sigma_beta;
}
transformed parameters {
  real<lower=0.0,upper=1.0> S[N_obs];
  real<lower=0.0,upper=1.0> E[N_obs];
  real<lower=0.0,upper=1.0> I[N_obs];
  real<lower=0.0,upper=1.0> R[N_obs];
  real<lower=0.0,upper=1.0> ptest[N_obs];
  //real R0[N_obs];
  E[1] = E0;
  S[1] = 0.0;
  I[1] = 0.0;
  R[1] = 0.0;
  ptest[1] = 0.0;
  for(j in 2:N_obs){
	  // S[j] = fmax(S[j-1] - beta[j]*S[j-1]*I[j-1],1E-15);
	  // E[j] = fmax(E[j-1] + beta[j]*S[j-1]*I[j-1] - sigma*E[j-1],1E-15);
	  // I[j] = fmax(I[j-1] + sigma*E[j-1]  - gamma*I[j-1],1E-15);
	  // R[j] = fmax(R[j-1] + gamma*I[j-1],1E-15);
	  S[j] = fmax(S[j-1] - beta[j]*S[j-1]*I[j-1],1E-15);
	  E[j] = fmax(E[j-1] + beta[j]*S[j-1]*I[j-1] - (1.0/5.0)*E[j-1],1E-15);
	  I[j] = fmax(I[j-1] + (1.0/5.0)*E[j-1]  - (1.0/7.0)*I[j-1],1E-15);
	  R[j] = fmax(R[j-1] + (1.0/7.0)*I[j-1],1E-15);
	  // S[j] = fmax(S[j-1] - beta*S[j-1]*I[j-1],1E-15);
	  // E[j] = fmax(E[j-1] + beta*S[j-1]*I[j-1] - (1.0/5.0)*E[j-1],1E-15);
	  // I[j] = fmax(I[j-1] + (1.0/5.0)*E[j-1]  - (1.0/7.0)*I[j-1],1E-15);
	  // R[j] = fmax(R[j-1] + (1.0/7.0)*I[j-1],1E-15);
	  ptest[j] = fmax(I[j]/(I[j] + s*(S[j] + E[j] + R[j])),1E-15);
  }  
  // for(i in 1:N_obs){
	  // R0[i] = beta[i]/gamma;
  // }
}
model {
  //gamma ~ normal(1.0/7.0,0.001);
  //sigma ~ normal(1.0/5.0,0.001);
  //beta ~ normal(0.4,0.5);
  beta[1] ~ normal(0.4,0.5);
  for(i in 2:N_obs){
	  beta[i] ~ normal(0.95*beta[i-1],0.01*beta[i]);
  }
  for(j in 2:N_obs){
	y[j] ~ poisson(fmax(Ntest[j]*ptest[j],1E-15));
  }
}
//generated quantities {
  // real x_pred[N_obs,4]; 
  // real y_pred[N_obs,4]; 
  // x_pred = integrate_ode_rk45(SIER, x0, 0, t_obs, theta, rep_array(0.0,0), rep_array(0,0),1e-8,1e-8,1e6) ;
  // for(i in 1:4){
	  // for(j in 1:N_obs){
		  // y_pred[j,i] = poisson_rng(fmax(x[j,i]*POP,0.0));
	  // }
  // }
// }


