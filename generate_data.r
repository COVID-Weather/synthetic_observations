library(rstan)
options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

setwd('d:/dropbox/working/covid19/synthetic_observations/' )

niter <- 100
dt    <- 1
S=E=I=R=t <- numeric(niter/dt)
I[1] <- 0.001
S[1] <- 1.0 - I[1]

gamma <- 0.1
beta <- 0.6
sigma <- 0.1

for(i in 2:niter){
	dSdt <- -beta*S[i-1]*I[i-1]
	dEdt <- beta*S[i-1]*I[i-1] - sigma*E[i-1]
	dIdt <- sigma*E[i-1] - gamma*I[i-1]
	dRdt <- gamma*I[i-1]
	
	S[i] <- S[i-1] + dSdt*dt
	E[i] <- E[i-1] + dEdt*dt
	I[i] <- I[i-1] + dIdt*dt
	R[i] <- R[i-1] + dRdt*dt
}

##--MAKE A PLOT OF MODEL SOLUTION--#######
plot(S,type='l',ylim=c(0,1))
lines(E)
lines(I)
lines(R)

##--GENERATE OBSERVATIONS-############
POP <- 1000              # total population size
t_obs <- seq(1,niter,1)  # observation points - once per day
N_obs <- length(t_obs)   # total number of obs

I_obs0 <- rpois(N_obs,lambda=I[t_obs]*POP)
S_obs0 <- rpois(N_obs,lambda=S[t_obs]*POP)
E_obs0 <- rpois(N_obs,lambda=E[t_obs]*POP)
R_obs0 <- rpois(N_obs,lambda=R[t_obs]*POP)

##--MAKE A PLOT OF OBSERVATIONS--##############
par(mfrow=c(2,2),mar=c(2,2,2,2))
plot(S,type='l')
	par(new=TRUE)
	plot(S_obs0,yaxt='n',xaxt='n')
	axis(side=4)
plot(E,type='l')
	par(new=TRUE)
	plot(E_obs0,yaxt='n',xaxt='n')
	axis(side=4)
plot(I,type='l')
	par(new=TRUE)
	plot(I_obs0,yaxt='n',xaxt='n')
	axis(side=4)
plot(R,type='l')
	par(new=TRUE)
	plot(R_obs0,yaxt='n',xaxt='n')
	axis(side=4)

#################################################################
## ESTIMATE PARAMETERS ##########################################
#################################################################
I0 <- 0.001
S0 <- 1.0 - I0
R0 <- 0.0
E0 <- 0.0

x0 <- c(S0,E0,I0,R0) #initial conditions

data <- list(POP=POP,y=cbind(S_obs0,E_obs0,I_obs0,R_obs0),N_obs=N_obs,t_obs=t_obs,x0=x0)

mod_SEIR <- stan_model('stan_SIER.stan') #compile the stan code
  
mcmc_SEIR <- sampling(mod_SEIR,data=data,open_progress=TRUE) #mcmc

post_SEIR <- extract(mcmc_SEIR) #extract samples

##--PLOT THE FIT--###################
par(mfrow=c(1,1))
y_pred <- apply(post_SEIR$y_pred,c(2,3),mean) 
matplot(y_pred,type='l',lty=1,ylim=c(0,1200))
for(i in 1:4){
	lines(apply(post_SEIR$y_pred[,,i],2,function(x) quantile(x,p=0.025)),col=i,lty=2)
	lines(apply(post_SEIR$y_pred[,,i],2,function(x) quantile(x,p=0.975)),col=i,lty=2)
	points(data$y[,i],col=i)
}
legend('top',legend=c('S','E','I','R'),col=1:4,lty=1,bty='n')



