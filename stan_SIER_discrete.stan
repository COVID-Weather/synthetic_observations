data {
  real<lower=0> POP;           //total population size		
  int<lower=0>  N_obs;         // num obs
  int           x_p;         // number of state variables
  int           N_obsvar;    // number of observed variables
  int<lower=0>  y[N_obs,x_p]; // observed variable at measurement times
  int           theta_p;  //number of fitted parameters
  int           i_obsvar[N_obsvar]; //indices for observations
  //int           i_obsvar; //indices for observations
  real          x0[x_p];
//  int           N_pred;
//  real          t_pred[N_pred];  
}
parameters {
  real<lower=0.0> beta;
  real<lower=0.0> gamma;
  real<lower=0.0> sigma;
  real x[N_obs,x_p];
  //real<lower=0.0,upper=1.0> s0;
}
model {
  beta ~ normal(0.5,2);
  gamma ~ normal(0.5,2);
  sigma ~ normal(0.5,2);
  sigma_mod ~ uniform(1E-15,1);
  
  x[1,] ~ poisson(x0*POP);
  for(i in 2:N_obs){
	  x[i,1] ~ poisson(x[i,1] - beta*x[i,1]*x[i,3], sigma_mod);
	  x[i,2] ~ normal(x[i,2] + beta*x[i,1]*x[i,3] - sigma*x[i,2], sigma_mod);
	  x[i,3] ~ normal(x[i,3] + sigma*x[i,2] - gamma*x[i,3], sigma_mod);
	  x[i,4] ~ normal(x[i,4] + gamma*x[i,3], sigma_mod);
  }
  
  for(i in i_obsvar){
	  for(j in 1:N_obs){
		y[j,i] ~ poisson(fmax(x[j,i]*POP,1E-15));
	  }
  }
}
// generated quantities {
  // real x_pred[N_pred,x_p]; 
  // real y_pred[N_pred,x_p];
  // x_pred = integrate_ode_rk45(SIER, x0, 0, t_pred, theta, rep_array(0.0,0), rep_array(0,0),1e-8,1e-8,1e6);
  // for(i in 1:x_p){
	  // for(j in 1:N_obs){
		  // y_pred[j,i] = poisson_rng(fmax(x[j,i]*POP,1E-15));
	  // }
  // }
// }


