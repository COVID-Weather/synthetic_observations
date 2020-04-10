library(rstan)
options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

setwd('d:/dropbox/working/covid19/synthetic_observations/' )

##############################################################
## GENERATE DATA #############################################
##############################################################
niter <- 100
dt    <- 1
S=E=I=R=t <- numeric(niter/dt)
I[1] <- 0.001
S[1] <- 1.0 - I[1]

gamma <- 0.2
beta <- 0.6
sigma <- 0.2

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
plot(S,type='l',ylim=c(0,1)); 
	par(new=TRUE); 
	#plot(S_obs0,yaxt='n',xaxt='n',ylim=c(0,990))
	#axis(side=4)
	mtext('Susceptible')
plot(E,type='l')
	#par(new=TRUE)
	#plot(E_obs0,yaxt='n',xaxt='n')
	#axis(side=4)
	mtext('Exposed')
plot(I,type='l')
	par(new=TRUE)
	plot(I_obs0,yaxt='n',xaxt='n')
	axis(side=4)
	mtext('Infected')
plot(R,type='l')
	par(new=TRUE)
	plot(R_obs0,yaxt='n',xaxt='n',ylim=c(0,990))
	axis(side=4)
	mtext('Recovered')

#################################################################
## ESTIMATE PARAMETERS ##########################################
#################################################################
I0 <- 0.001
S0 <- 1.0 - I0
R0 <- 0.0
E0 <- 0.0
x0 <- c(S0,E0,I0,R0) #initial conditions

num <- 0.6
data <- list(POP=POP,
		     y=cbind(S_obs0,E_obs0,as.integer(num*I_obs0),as.integer(num*R_obs0)),
			 N_obs=N_obs,
			 t_obs=t_obs,
			 x0 = x0,
			 x_p=4,
			 N_obsvar=1,
			 theta_p=3,
			 i_obsvar=c(3))
opt_SEIR <- optimizing(mod_SEIR,data=data)
opt_SEIR$par[1:3]
opt_SEIR$par[2]/opt_SEIR$par[1]




##--COMPILE STAN CODE--###############
mod_SEIR <- stan_model('stan_SIER.stan') #compile the stan code
#mod_SEIR <- stan_model('stan_SIER_multinomial.stan')
  
##--MCMC VS OPTIMIZATION--############
mcmc_SEIR <- sampling(mod_SEIR,data=data,open_progress=TRUE) #mcmc

post_SEIR <- extract(mcmc_SEIR) #extract samples
post_SEIR_I_R <- extract(mcmc_SEIR_I_R)


####################################################################
## LINEAR BIAS wrt TIME ############################################
## - code assumes a linear ramp-up in the proportion of infected being tested
## - p0 is the initial proportion of infecteds being observed 
####################################################################
data_I_R <- list(POP=POP,
		     y=cbind(S_obs0,E_obs0,I_obs0,R_obs0),
			 N_obs=N_obs,
			 t_obs=t_obs,
			 x0 = x0,
			 x_p=4,
			 N_obsvar=2,
			 theta_p=3,
			 i_obsvar=c(3,4))

opt_SEIR <- optimizing(mod_SEIR,data=data_I_R)

n_sim <- 20
p0s <- seq(0,0.9,length.out=n_sim)
PARS <- matrix(NA,nrow=3,ncol=n_sim)
SDS  <- matrix(NA,nrow=3,ncol=n_sim)

for(i in 1:20){
	data_I_R_test <- data_I_R
	data_I_R_test$y[,3] <- as.integer(seq(p0s[i],1,length.out=N_obs)*data_I_R_test$y[,3])
	
	opt_SEIR_I_R_test <- optimizing(mod_SEIR,data=data_I_R_test,init='0',hessian=TRUE)
	PARS[,i] <- opt_SEIR_I_R_test$par[1:3]
	SDS[,i] <- sqrt(diag(solve(-opt_SEIR_I_R_test$hessian)))
}


####################################################################
## PIECE-WISE BIAS #################################################
####################################################################
n_sim_pw <- 20
p0s_pw <- seq(0,0.9,length.out=n_sim_pw)
PARS_pw <- matrix(NA,nrow=3,ncol=n_sim_pw)
SDS_pw  <- matrix(NA,nrow=3,ncol=n_sim_pw)

for(i in 1:20){
	data_I_R_test <- data_I_R
	data_I_R_test$y[,3] <- as.integer(c(rep(p0s[i],10),rep(1,10))*data_I_R_test$y[,3])
	
	opt_SEIR_I_R_test <- optimizing(mod_SEIR,data=data_I_R_test,init='0',hessian=TRUE)
	PARS_pw[,i] <- opt_SEIR_I_R_test$par[1:3]
	SDS_pw[,i] <- sqrt(diag(solve(-opt_SEIR_I_R_test$hessian)))
}


par(mfrow=c(2,2),mar=c(3,3,3,3),oma=c(2,2,2,2))
plot(seq(p0s[1],1,length.out=N_obs),type='l')
for(i in 2:length(p0s)) lines(seq(p0s[i],1,length.out=N_obs))
mtext(side=2,'Testing bias',line=2.5)
plot(p0s,PARS[1,],ylim=c(0.1,1.1),type='l',col='blue')
	lines(p0s,PARS[1,] + 2*SDS[1,],lty=2,col='blue')
	lines(p0s,PARS[1,] - 2*SDS[1,],lty=2,col='blue')
	legend('topright',legend=c(expression(italic(gamma)),
	                           expression(italic(beta)),
							   expression(italic(sigma)),
							   expression(italic(R['0']))),bty='n',lty=1,col=c('blue','red','orange','black'))
lines(p0s,PARS[2,],col='red')
	lines(p0s,PARS[2,] + 2*SDS[2,],lty=2,col='red')
	lines(p0s,PARS[2,] - 2*SDS[2,],lty=2,col='red')
abline(h=gamma,col='blue')
abline(h=beta,col='red')
lines(p0s,PARS[3,],col='orange')
	lines(p0s,PARS[3,] + 2*SDS[3,],lty=2,col='orange')
	lines(p0s,PARS[3,] - 2*SDS[3,],lty=2,col='orange')
abline(h=sigma-0.05,col='orange')

par(new=TRUE)
# plot(p0s,PARS[2,]/PARS[1,],ylim=c(2.5,4),xaxt='n',yaxt='n',type='l')
	# lines(p0s,PARS[2,]/PARS[1,] + 2*sqrt(SDS[1,]^2 + SDS[2,]^2),lty=2)
	# lines(p0s,PARS[2,]/PARS[1,] - 2*sqrt(SDS[1,]^2 + SDS[2,]^2),lty=2)
plot(p0s_pw,PARS_pw[2,]/PARS_pw[1,],ylim=c(2.5,4),xaxt='n',yaxt='n',type='l')
	lines(p0s,PARS_pw[2,]/PARS_pw[1,] + 2*sqrt(SDS_pw[1,]^2 + SDS_pw[2,]^2),lty=2)
	lines(p0s,PARS_pw[2,]/PARS_pw[1,] - 2*sqrt(SDS_pw[1,]^2 + SDS_pw[2,]^2),lty=2)
axis(side=4)
abline(h=beta/gamma)
mtext(side=1,'Initial testing bias',line=2.5)



#par(mfrow=c(1,2),mar=c(3,3,3,3),oma=c(2,2,2,2))
plot(c(rep(p0s_pw[1],50),rep(1,50)),type='l')
for(i in 2:length(p0s_pw)) lines(c(rep(p0s_pw[i],50),rep(1,50)))
mtext(side=2,'Testing bias',line=2.5)

plot(p0s_pw,PARS_pw[1,],ylim=c(0.1,1.1),type='l',col='blue')
	lines(p0s_pw,PARS_pw[1,] + 2*SDS_pw[1,],lty=2,col='blue')
	lines(p0s_pw,PARS_pw[1,] - 2*SDS_pw[1,],lty=2,col='blue')
	legend('topright',legend=c(expression(italic(gamma)),
	                           expression(italic(beta)),
							   expression(italic(sigma)),
							   expression(italic(R['0']))),bty='n',lty=1,col=c('blue','red','orange','black'))
lines(p0s_pw,PARS_pw[2,],col='red')
	lines(p0s_pw,PARS_pw[2,] + 2*SDS_pw[2,],lty=2,col='red')
	lines(p0s_pw,PARS_pw[2,] - 2*SDS_pw[2,],lty=2,col='red')
lines(p0s,PARS_pw[3,],col='orange')
	lines(p0s,PARS_pw[3,] + 2*SDS_pw[3,],lty=2,col='orange')
	lines(p0s,PARS_pw[3,] - 2*SDS_pw[3,],lty=2,col='orange')
abline(h=sigma-0.05,col='orange')

abline(h=gamma,col='blue')
abline(h=beta,col='red')
par(new=TRUE)
plot(p0s_pw,PARS_pw[2,]/PARS_pw[1,],ylim=c(2.5,4),xaxt='n',yaxt='n',type='l')
	lines(p0s,PARS_pw[2,]/PARS_pw[1,] + 2*sqrt(SDS_pw[1,]^2 + SDS_pw[2,]^2),lty=2)
	lines(p0s,PARS_pw[2,]/PARS_pw[1,] - 2*sqrt(SDS_pw[1,]^2 + SDS_pw[2,]^2),lty=2)
axis(side=4)
abline(h=beta/gamma)
mtext(side=1,'Initial testing bias',line=2.5)

