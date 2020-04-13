data {
  real<lower=0> POP;          	
  int<lower=0>  N_obs;         
  int<lower=0>  y[N_obs]; 
  real          x0[4];
}
parameters {
  real<lower=0.0> theta[1];
}
transformed parameters {
  real x[N_obs,4];
  x[1,] = x0;
  for(j in 2:N_obs){
	  x[j,1] = x[j-1,1] - theta[1]*x[j-1,1]*x[j-1,3];
	  x[j,2] = x[j-1,2] + theta[1]*x[j-1,1]*x[j-1,3] - (1.0/5.0)*x[j-1,2];
	  x[j,3] = x[j-1,3] + (1.0/5.0)*x[j-1,2]  - (1.0/7.0)*x[j-1,3];
	  x[j,4] = x[j-1,4] + (1.0/7.0)*x[j-1,3];
  }  
}
model {
  theta[1] ~ normal(0.2,2);
	  for(j in 2:N_obs){
		y[j] ~ poisson((x[j,3]-x[j-1,3])*POP);
	  }
}
// generated quantities {
  // real x_pred[N_obs,4]; 
  // real y_pred[N_obs,4];
  // x_pred = integrate_ode_rk45(SIER, x0, 0, t_obs, theta, rep_array(0.0,0), rep_array(0,0),1e-8,1e-8,1e6) ;
  // for(i in 1:4){
	  // for(j in 1:N_obs){
		  // y_pred[j,i] = poisson_rng(fmax(x[j,i]*POP,0.0));
	  // }
  // }
// }


