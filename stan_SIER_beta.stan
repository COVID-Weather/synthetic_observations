functions {
  real[] SIER(real   t,           // time
               real[] x,           // state x[1]:CH  x[2]:PR, x[3]:Chl , x[4]:N
               real[] theta,
               real[] x_r,
               int[]  x_i) {       // parameters

    real beta  = theta[1]; 
    
    real S = fmax(x[1],0.0);
    real E = fmax(x[2],0.0);
    real I = fmax(x[3],0.0);
    real R = fmax(x[4],0.0);

	real dS = -beta*S*I;
	real dE =  beta*S*I - (1.0/14.0)*E;
	real dI = (1.0/14.0)*E - (1.0/10.0)*I;
	real dR = (1.0/10.0)*I;

    return {dS,dE,dI,dR};
  }
}
data {
  real<lower=0> POP;           //total population size		
  int<lower=0>  N_obs;         // num obs
  real          t_obs[N_obs];          // obs times
  int<lower=0>  y[N_obs,4]; // observed variable at measurement times
  real          x0[4];
}
parameters {
  real<lower=0.0> theta[1];
}
transformed parameters {	
  real x[N_obs,4] = integrate_ode_rk45(SIER, x0, 0, t_obs, theta, rep_array(0.0,0), rep_array(0,0), 1e-8, 1e-8, 1e6);
}
model {
  theta[1] ~ normal(0.2,1);
  for(i in 1:4){
	  //y[1:N_obs,i] ~ poisson(x[1:N_obs,i]*POP)
	  for(j in 1:N_obs){
		y[j,i] ~ poisson(fmax(x[j,i]*POP,0.0));
	  }
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


