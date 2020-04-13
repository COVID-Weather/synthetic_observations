library(rstan)
options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
library(colorRamps)

setwd('d:/dropbox/working/covid19/synthetic_observations/' )

##############################################################
## READ DATA #################################################
##############################################################
dat    <- read.csv('daily.csv',stringsAsFactor=FALSE)
abbrev <- read.csv('state_abbrev.csv',stringsAsFactor=FALSE)
pops   <- read.csv('co-est2019-alldata.csv',stringsAsFactor=FALSE)

states <- unique(dat$state)

E0 <- 0.00002
S0 <- 1.0 - E0
R0 <- 0.0
I0 <- 0.0

x0 <- c(S0,E0,I0,R0) #initial conditions

########################################################
## TESTING EXPERIMENTS #################################
########################################################
fs <- seq(0.25,2.5,length.out=10)
cols <- terrain.colors(20)

Nout <- 65:105
par(mfrow=c(1,1),mar=c(2,2,2,2),oma=c(2,2,2,2))
plot(-999,xlim=range(Nout),ylim=c(0,3),bty='n')
for(i in 1:length(fs)){
	p <- c(rep(1,0),seq(1,fs[i],length.out=length(Nout)-0))
	lines(Nout,p,lty=2,col=cols[i])
	abline(h=1)
}
mtext(side=2,line=2.5,'Case Multiplier')

#############################
## PLOT #####################
#############################
par(mfrow=c(8,7),mar=c(0,2,0,0),oma=c(4,2,1,1) ,cex.axis=0.7)
DATA <- list()
k <- 1
#for(i in c(1:3,5:56)){
for(i in 1:length(states)){

	states[i]
	name <- abbrev$name[abbrev$abbrev==states[i]]	
	POP <- pops$POPESTIMATE2019[pops$STNAME==name & pops$COUNTY==0]
	
	st <- dat[dat$state==states[i],]
	st$doy <- as.numeric(strftime(paste(substr(as.character(st$date),1,4), 
										 substr(as.character(st$date),5,6), 
										 substr(as.character(st$date),7,8), sep='-'), format = "%j")) 
	dPOS <- rev(st$positiveIncrease);     dPOS[1] <- 0
	dNEG <- rev(st$negativeIncrease);     dNEG[1] <- 0
	doy  <- rev(st$doy)
	dHOS <- rev(st$hospitalizedIncrease); dHOS[1] <- 0

	N_obs <- length(dPOS)
	
	d <- dPOS 
	d[d<0] <- 0
	if(length(POP)>0){
	DATA[[k]] <- list(POP=POP,
				      y=d,
					  dPOS=dPOS,
					  dNEG=dNEG,
				      N_obs=N_obs, 
				      x0 = x0,
				      doy=doy,
					  state=states[i])
	k <- k+1				  
	}

	plot(doy,dPOS,type='l',xlim=c(65,102),xaxt='n',bty='n')
	mtext(states[i],line=-1.5,adj=0.05)
	if(i>49) axis(side=1)
}
mtext(side=2,outer=TRUE,'Daily Cases')
mtext(side=1,outer=TRUE,'Day of Year',line=2.5)



par(mfrow=c(8,7),mar=c(0,2,0,0),oma=c(3,2,1,1),cex.axis=0.7)
for(i in 1:length(DATA)){
	plot(DATA[[i]]$doy,DATA[[i]]$dPOS + DATA[[i]]$dNEG,xaxt='n',type='l',bty='n',xlim=c(63,102))
	lw <- loess(DATA[[i]]$dPOS + DATA[[i]]$dNEG  ~ DATA[[i]]$doy)
	lines(DATA[[i]]$doy,lw$fitted,lwd=1.5,col='red')
	legend('topleft',legend=states[i],bty='n')
	if(i>44) axis(side=1)
}
mtext(side=1,outer=TRUE,'Day of Year',line=1)
mtext(side=2,outer=TRUE,'Number of Tests',line=0.5)

par(mfrow=c(8,7),mar=c(0,0,0,0),oma=c(3,2,1,1),cex.axis=0.7)
for(i in 1:length(DATA)){
	dat <- DATA[[i]]
	plot(dat$dPOS+dat$dNEG,dat$dPOS,xaxt='n',yaxt='n',ylim=c(),xlim=c(),pch=19)
	legend('topleft',legend=states[i],bty='n',cex=1.2)
}
mtext(side=1,outer=TRUE,'Number of Tests')
mtext(side=2,outer=TRUE,'Number of Positives',line=0.5)
##################################################################
## FIT MODELS ####################################################
##################################################################
############################
##--COMPILE--###############
############################
#mod_b   <- stan_model('stan_SIER_beta.stan') #compile the stan code
#mod_b_s <- stan_model('stan_SIER_beta_sigma.stan') #compile the stan code
#mod     <- stan_model('stan_SEIR_discrete.stan')
mod_betaAR1 <- stan_model('stan_SEIR_discrete_beta_AR1.stan')

##--TRY IT OUT--####
# i=2
# opt <- optimizing(mod_betaAR1,data=DATA[[i]],algorithm="Newton",hessian=TRUE,init='0')
# mcmc <- sampling(mod_betaAR1,data=DATA[[i]],open_progress=TRUE,init='0')
# post <- extract(mcmc)

# q <- apply(post$R0,2,function(x) quantile(x, p=c(0.05,0.5,0.95)))
# matplot(DATA[[i]]$doy,t(q),type='l',lty=c(2,1,2),col='black')

#############################
## LOOP OVER STATES #########
#############################


#ex <- c(4,20,28,53,56)
#ex <- c(4,7,13,25,28,30,53,56)
ex <- c(4,23)

FITS <- list()
par(mfrow=c(8,7),mar=c(0,2,0,0),oma=c(1,1,1,1))
for(i in c(1:56)[-ex]){
	opt <-try(optimizing(mod_betaAR1,data=DATA[[i]],hessian=TRUE,algorithm='Newton',init='0'),silent=TRUE)
	if(class(opt)!='try-error'){
	FITS[[DATA[[i]]$state]] <- opt
		gamma <- opt$par[names(opt$par)=='gamma']
		N <- length(DATA[[i]]$y)
		plot(DATA[[i]]$doy,opt$par[1:N]/gamma,xaxt='n',ylim=c(0,10),type='l',bty='n',xlim=c(63,102))
		abline(h=1,lty=2)
		#plot(opt$par[1:38])
		sd_gamma <- sqrt(diag(solve(-opt$hessian)))[N+1]
		var <- try((diag(solve(-opt$hessian))[1:N]) + sd_gamma^2)
		lines(DATA[[i]]$doy,opt$par[1:N]/gamma + 4*sqrt(var))
		lines(DATA[[i]]$doy,opt$par[1:N]/gamma - 4*sqrt(var))
		mtext(DATA[[i]]$state,line=-1.5,adj=0.1)
	}
}

# MCMC <- list()
# par(mfrow=c(8,7),mar=c(0,2,0,0),oma=c(1,1,1,1))
# for(i in c(1:56)[-ex]){
# print(i)
	# mcmc <- sampling(mod_betaAR1,data=DATA[[i]],init='0')
	# post <- extract(mcmc)
	# MCMC[[DATA[[i]]$state]] <- mcmc
	
# }


##--MULTIPLE LINES--############
fs <- seq(0.2,2.5,length.out=5)
cols <- terrain.colors(7)

par(mfrow=c(8,7),mar=c(0,2,0,0),oma=c(1,1,1,1))
for(i in c(1:56)[-ex]){
	
	data <- DATA[[i]]
	dPOS <- data$y
	N_obs <- data$N_obs
	
	opt <- try(optimizing(mod_betaAR1,data=data,algorithm="Newton",hessian=TRUE),silent=TRUE)
	if(class(opt)!='try-error'){
		gamma <- opt$par[names(opt$par)=='gamma']
		N <- length(DATA[[i]]$y)
		plot(DATA[[i]]$doy,opt$par[1:N]/gamma,xaxt='n',ylim=c(0,10),type='l',bty='n',xlim=c(63,102))
		abline(h=1,lty=2)
		
		#sd_gamma <- sqrt(diag(solve(-opt$hessian)))[N+1]
		#var <- (diag(solve(-opt$hessian))[1:N]) + sd_gamma^2
		#lines(DATA[[i]]$doy,opt$par[1:N]/gamma + 4*sqrt(var))
		#lines(DATA[[i]]$doy,opt$par[1:N]/gamma - 4*sqrt(var))
		mtext(DATA[[i]]$state,line=-1.5,adj=0.1)
		
		for(j in 1:length(fs)){
			p <- seq(1,fs[j],length.out=N_obs)	
			data$y <- as.integer(apply(matrix(dPOS*p),1,function(x) max(x,0)))
		
			opt <- try(optimizing(mod_betaAR1,data=data,algorithm="Newton",hessian=TRUE,init='0'),silent=TRUE)
			if(class(opt)!='try-error'){
				gamma <- opt$par[names(opt$par)=='gamma']
				lines(data$doy,opt$par[1:N]/gamma,col=cols[j])
			}
		}
}
}


