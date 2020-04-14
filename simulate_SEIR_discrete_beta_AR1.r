
niter <- 200
nsamp <- 200

Npop <- 1E4
I=S=E=R <- matrix(0,ncol=niter,nrow=nsamp)
I0=S0=E0=R0 <- matrix(0,ncol=niter,nrow=nsamp)

Sobs=Eobs=Iobs=Robs <- matrix(0,ncol=niter,nrow=nsamp)
S0obs=E0obs=I0obs=R0obs <- matrix(0,ncol=niter,nrow=nsamp)

E[,1]=E0[,1] <- 0.001
S[,1]=S0[,1] <- 1 - E[,1]

gamma <- 0.2
sigma <- 0.1
beta=beta0 <- matrix(NA,ncol=niter,nrow=nsamp)

rw_sd <- 0.02
rw_slope <- -0.005
beta_init <- 0.7
beta_init_sd <- 0.2

beta[,1] <- rnorm(nsamp,mean=beta_init,sd=beta_init_sd)
beta0[,1] <- rnorm(nsamp,mean=beta_init,sd=beta_init_sd)
#beta[,1] <- rnorm(nsamp,mean=-0.2,sd=0.3)  #exponent
#beta0[,1] <- rnorm(nsamp,mean=-0.2,sd=0.3)

for(i in 2:niter){
	# S[,i] = S[,i-1] - beta[,i-1]*I[,i-1]*S[,i-1]
	# E[,i] = E[,i-1] + beta[,i-1]*I[,i-1]*S[,i-1] - sigma*E[,i-1]
	S[,i] = S[,i-1] - beta[,i-1]*I[,i-1]*S[,i-1]
	E[,i] = E[,i-1] + beta[,i-1]*I[,i-1]*S[,i-1] - sigma*E[,i-1]
	I[,i] = I[,i-1] + sigma*E[,i-1] - gamma*I[,i-1]
	R[,i] = R[,i-1] + gamma*I[,i-1]
	
	# Sobs[,i] <- rpois(nsamp,lambda=S[,i]*Npop)
	# Eobs[,i] <- rpois(nsamp,lambda=E[,i]*Npop)
	# Iobs[,i] <- rpois(nsamp,lambda=I[,i]*Npop)
	# Robs[,i] <- rpois(nsamp,lambda=R[,i]*Npop)
		
	#beta[,i] <- beta[,i-1] -0.015 + rnorm(nsamp, mean=0, sd=0.05)
	beta[,i] = apply(matrix(beta[,i-1] + rnorm(nsamp,mean=rw_slope,sd=rw_sd)), 1, function(x) max(0,x))

	S0[,i] = S0[,i-1] - beta0[,i-1]*I0[,i-1]*S0[,i-1]
	E0[,i] = E0[,i-1] + beta0[,i-1]*I0[,i-1]*S0[,i-1] - sigma*E0[,i-1]
	I0[,i] = I0[,i-1] + sigma*E0[,i-1] - gamma*I0[,i-1]
	R0[,i] = R0[,i-1] + gamma*I0[,i-1]

	# S0obs[,i] <- rpois(nsamp,lambda=S0[,i]*Npop)
	# E0obs[,i] <- rpois(nsamp,lambda=E0[,i]*Npop)
	# I0obs[,i] <- rpois(nsamp,lambda=I0[,i]*Npop)
	# R0obs[,i] <- rpois(nsamp,lambda=R0[,i]*Npop)
	# beta0[,i] <- beta0[,i-1] + rnorm(nsamp, mean=0, sd=0.02)
	beta0[,i] = apply(matrix(beta0[,i-1] + rnorm(nsamp,mean=0,sd=rw_sd)), 1, function(x) max(0,x))
}

#beta <- exp(beta)
#beta0 <- exp(beta0)

#####################################################
## PLOT SIMULATION ##################################
#####################################################
##-TIME-VARYING BETA-##############
#pdf('')
par(mfrow=c(4,3),mar=c(2,2,2,2))
hist(beta[,1]/(1/7),breaks=20,main='')
mtext(expression('R'['0']*'(t'['0']*')'))

matplot(t(beta)/(1/7),type='l',bty='n')
	mtext(side=2,expression('R'['0']),line=2.5)
	#mtext(expression(italic(beta['t']*' = f'['cont']*'(t) * p(E|cont.)')))
plot(apply(beta/(1/7),2,function(x) quantile(x,p=0.05)),type='l',ylim=c(0,10),ylab='t(beta)')
	lines(apply(beta/(1/7),2,function(x) quantile(x,p=0.5)))
	lines(apply(beta/(1/7),2,function(x) quantile(x,p=0.95)))
	mtext(expression(italic('p'['0.05, 0.5, 0.95']*'('*beta['t']*')')))

matplot(t(S),type='l')
	mtext('Susceptible')
#matplot(t(E),type='l')
#	mtext('Exposed')
matplot(t(I),type='l')
	mtext('Infected')
matplot(t(R),type='l')
	mtext('Recovered')

hist(beta0[,1]/(1/7),breaks=20,main='')
mtext(expression('R'['0']*'(t'['0']*')'))

matplot(t(beta0)/(1/7),type='l',bty='n')
	mtext(side=2,expression('R'['0']),line=2.5)
	#mtext(expression(italic(beta['t']*' = f'['cont']*'(t) * p(E|cont.)')))
plot(apply(beta0/(1/7),2,function(x) quantile(x,p=0.05)),type='l',ylim=c(0,10),ylab='t(beta)')
	lines(apply(beta/(1/7),2,function(x) quantile(x,p=0.5)))
	lines(apply(beta/(1/7),2,function(x) quantile(x,p=0.95)))
	mtext(expression(italic('p'['0.05, 0.5, 0.95']*'('*beta['t']*')')))


matplot(t(S0),type='l')
	mtext('Susceptible')
#matplot(t(E),type='l')
#	mtext('Exposed')
matplot(t(I0),type='l')
	mtext('Infected')
matplot(t(R0),type='l')
	mtext('Recovered')


##--CUMULATIVE CASES--####################
par(mfrow=c(1,2),mar=c(3,3,3,3),oma=c(2,2,2,2))
plot(-1E5,xlim=c(1,70),ylim=c(100,5E4),log='y')
for(i in 1:100){ 
	tmp <- cumsum(Iobs[i,])
	tmpp <- tmp[tmp>100]
	lines(tmpp,col=i)
}
lines(100*exp(0.01*seq(1,100)),lty=2,lwd=3)
lines(100*exp(0.1*seq(0,100)),lty=2,lwd=3)
lines(100*exp(1*seq(0,100)),lty=2,lwd=3)
mtext(side=1,'Days since 100th case',line=2.5)	
mtext(side=2,'Cumulative cases',line=2.5)	

plot(-1E5,xlim=c(1,150),ylim=c(100,5E4),log='y')
for(i in 1:100){ 
	tmp <- cumsum(I0obs[i,])
	tmpp <- tmp[tmp>100]
	lines(tmpp,col=i)
}
lines(100*exp(0.01*seq(1,150)),lty=2,lwd=3)
lines(100*exp(0.1*seq(0,150)),lty=2,lwd=3)
lines(100*exp(1*seq(0,150)),lty=2,lwd=3)
mtext(side=1,'Days since 100th case',line=2.5)	
mtext(side=2,'Cumulative cases',line=2.5)	



