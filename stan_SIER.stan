functions {
  real[] SIER(real   t,           // time
               real[] x,           // state x[1]:CH  x[2]:PR, x[3]:Chl , x[4]:N
               real[] theta,
               real[] x_r,
               int[]  x_i) {       // parameters

	
    //real gamma = theta[1];    
    real beta  = theta[1]; 
    real sigma = theta[2]; 
    
    real S = fmax(x[1],0.0);
    real E = fmax(x[2],0.0);
    real I = fmax(x[3],0.0);
    real R = fmax(x[4],0.0);

	real dS = -beta*S*I;
	real dE =  beta*S*I - sigma*E;
	real dI = sigma*E - 0.1*I;
	real dR = 0.1*I;

    return {dS,dE,dI,dR};
  }
}
data {
  real<lower=0> POP;           //total population size		
  int<lower=0>  N_obs;         // num obs
  real          t_obs[N_obs];          // obs times
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
  real<lower=0.0>           theta[2];   // parameters
  //real<lower=0.0,upper=1.0> s0;
}
transformed parameters {
  real x[N_obs,x_p];
  x = integrate_ode_rk45(SIER, x0, 0, t_obs, theta, 
		rep_array(0.0,0), rep_array(0,0), 1e-8, 1e-8, 1e6);
}
model {
  // s0       ~ uniform(0.5,1.0);
  //theta[1] ~ normal(0.3,0.5);
  theta[1] ~ normal(0.3,0.5);
  theta[2] ~ normal(0.3,0.5);
  
  //for(i in i_obsvar){
	  for(j in 1:N_obs){
		y[j,3] ~ poisson(fmax(x[j,3]*POP,1E-15));
	  }
  //}
}
generated quantities {
  real x_pred[N_pred,x_p]; 
  real y_pred[N_pred,x_p];
  x_pred = integrate_ode_rk45(SIER, x0, 0, t_pred, theta, rep_array(0.0,0), rep_array(0,0),1e-8,1e-8,1e6);
  for(i in 1:x_p){
	  for(j in 1:N_obs){
		  y_pred[j,i] = poisson_rng(fmax(x[j,i]*POP,1E-15));
	  }
  }
}


