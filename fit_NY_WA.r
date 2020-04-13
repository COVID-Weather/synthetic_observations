library(rstan)
options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
library(colorRamps)

setwd('d:/dropbox/working/covid19/synthetic_observations/' )

##############################################################
## READ DATA #################################################
##############################################################
dat <- read.csv('daily.csv')
NY <- dat[dat$state=='NY',]
NY$doy <- as.numeric(strftime(paste(substr(as.character(NY$date),1,4), 
							         substr(as.character(NY$date),5,6), 
									 substr(as.character(NY$date),7,8), sep='-'), format = "%j")) 
dPOSny <- rev(NY$positiveIncrease);     dPOSny[1] <- 0
dNEGny <- rev(NY$negativeIncrease);     dNEGny[1] <- 0
doyny  <- rev(NY$doy)
dHOSny <- rev(NY$hospitalizedIncrease); dHOSny[1] <- 0
POPny  <- 19.54E7

WA <- dat[dat$state=='WA',]
WA$doy <- as.numeric(strftime(paste(substr(as.character(WA$date),1,4), 
							         substr(as.character(WA$date),5,6), 
									 substr(as.character(WA$date),7,8), sep='-'), format = "%j")) 
dPOSwa <- rev(WA$positiveIncrease);     dPOSwa[1] <- 0
dNEGwa <- rev(WA$negativeIncrease);     dNEGwa[1] <- 0
doywa  <- rev(WA$doy)
dHOSwa <- rev(WA$hospitalizedIncrease); dHOSwa[1] <- 0
POPwa  <- 7.62E7

niterny <- length(dPOSny)
niterwa <- length(dPOSwa)
########################
## PLOT DATA ###########
########################
##-NY--###
par(mfrow=c(3,2),mar=c(2,2,2,2))
plot(doyny,dPOSny,type='l');                 mtext('Positive (#ind.)')
plot(doyny,dNEGny,type='l');                 mtext('Negative (#ind.)')
plot(doyny,dPOSny + dNEGny,type='l');          mtext('Total Tests (#)')
plot(doyny,dPOSny/(dPOSny+dNEGny),type='l');     mtext('Positive Rate (#+/#)')
plot(doyny,(dPOSny + dNEGny)/POPny,type='l') ;  mtext('Tests as Fraction of Population')
plot(doyny,cumsum(dPOSny),log='y',type='l'); mtext('log Cumulative Positives')
#plot(doy,dHOS,type='l');                 mtext('New Hospitalized')
#plot(doy,cumsum(dHOS),log='y',type='l'); mtext('log Cumulative Hospitalized')

##-WA--###
par(mfrow=c(3,2),mar=c(2,2,2,2))
plot(doywa,dPOSwa,type='l');                 mtext('Positive (#ind.)')
plot(doywa,dNEGwa,type='l');                 mtext('Negative (#ind.)')
plot(doywa,dPOSwa + dNEGwa,type='l');          mtext('Total Tests (#)')
plot(doywa,dPOSwa/(dPOSwa+dNEGwa),type='l');     mtext('Positive Rate (#+/#)')
plot(doywa,(dPOSwa + dNEGwa)/POPwa,type='l') ;  mtext('Tests as Fraction of Population')
plot(doywa,cumsum(dPOSwa),log='y',type='l'); mtext('log Cumulative Positives')
#plot(doy,dHOS,type='l');                 mtext('New Hospitalized')
#plot(doy,cumsum(dHOS),log='y',type='l'); mtext('log Cumulative Hospitalized')

#########################################################
## PLOTS DATA WITH HIGH AND LOW POSITIVE RATES ##########
#########################################################
nyup <- 0.75
nydn <- 0.1

waup <- 0.3
wadn <- 0.02

par(mfrow=c(2,2))
plot(dPOSny,type='l',ylim=c(0,20000),bty='n')
lines(nyup*(dPOSny+dNEGny),lty=2)
lines(nydn*(dPOSny+dNEGny),lty=2)

plot(dPOSwa,type='l',ylim=c(0,3000),bty='n')
lines(waup*(dPOSwa+dNEGwa),lty=2)
lines(wadn*(dPOSwa+dNEGwa),lty=2)

plot(cumsum(dPOSny),type='l',ylim=c(10,1000000),log='y',bty='n')
lines(cumsum(nyup*(dPOSny+dNEGny)),lty=2)
lines(cumsum(nydn*(dPOSny+dNEGny)),lty=2)

plot(cumsum(dPOSwa),type='l',log='y',ylim=c(10,22000),bty='n')
lines(cumsum(waup*(dPOSwa+dNEGwa)),lty=2)
lines(cumsum(wadn*(dPOSwa+dNEGwa)),lty=2)

mtext(outer=TRUE,side=1,'Day',cex=1.2)

##########################################################
##--PLOTS OF TESTING AND POSITITIVE RATES OVER TIME--#####
##########################################################
cols <- matlab.like(length(doyny))

par(mfrow=c(2,3))
plot(-999,xlim=c(64,99),ylim=c(0,12000),bty='n')
x <- doyny; y <- dPOSny
segments(x[-length(x)],y[-length(y)],x[-1L],y[-1L],col=cols)
points(x,y,col='black',bg=cols,pch=21,lwd=1)

plot(-999,xlim=c(0,25000),ylim=c(0,12000),bty='n')
x <- dPOSny + dNEGny; y <- dPOSny
segments(x[-length(x)],y[-length(y)],x[-1L],y[-1L],col=cols)
points(x,y,col='black',bg=cols,pch=21,lwd=1)

plot(-999,xlim=c(0,25000),ylim=c(0,1),bty='n')
x <- dPOSny + dNEGny; y <- dPOSny/(dPOSny + dNEGny)
y[y=='NaN'] <- 0; x[x=='NaN'] <- 0
segments(x[-length(x)],y[-length(y)],x[-1L],y[-1L],col=cols)
points(x,y,col='black',bg=cols,pch=21,lwd=1)

##--WASHINGTON--##########
plot(-999,xlim=c(64,99),ylim=c(0,1000),bty='n')
x <- doywa; y <- dPOSwa
segments(x[-length(x)],y[-length(y)],x[-1L],y[-1L],col=cols)
points(x,y,col='black',bg=cols,pch=21,lwd=1)

plot(-999,xlim=c(0,12000),ylim=c(0,1000),bty='n')
x <- dPOSwa + dNEGwa; y <- dPOSwa
y[y=='NaN'] <- 0; x[x=='NaN'] <- 0
segments(x[-length(x)],y[-length(y)],x[-1L],y[-1L],col=cols)
points(x,y,col='black',bg=cols,pch=21,lwd=1)

plot(-999,xlim=c(0,12000),ylim=c(0,1),bty='n')
x <- dPOSwa + dNEGwa; y <- dPOSwa/(dPOSwa + dNEGwa)
#y[y=='NaN'] <- 0; x[x=='NaN'] <- 0
segments(x[-length(x)],y[-length(y)],x[-1L],y[-1L],col=cols)
points(x,y,col='black',bg=cols,pch=21,lwd=1)
image.plot(legend.only=TRUE,matrix(doy))

#########################
## PRIOR ON R0 ##########
#########################
gamma <- 1/14
xin <- seq(0,15,0.01)
par(mfrow=c(1,1))
plot(xin,dnorm(xin,mean=0.2/gamma,sd=3))

###################################################
## PLOTS OF OBSERVING EXPERIMENTS #################
###################################################
par(mfrow=c(3,1),mar=c(2,2,2,2))
plot(-999,xlim=c(0,40),ylim=c(0,2.5),bty='n')
for(i in 1:length(fs)){
	p <- c(rep(1,0),seq(1,fs[i],length.out=niterny-0))
	
	lines(p,lty=2,col=cols[i])
}
	abline(h=1,lty=1,lwd=2)
	mtext(side=2,line=2.5,'Case Multiplier')
plot(doyny,dPOSny,ylim=c(0,3E4),type='l',lwd=2,bty='n',xlim=c(60,101))
fs <- seq(0.25,2.5,length.out=10)
yny <- matrix(NA,ncol=length(fs),nrow=niterny)
for(i in 1:length(fs)){
	p <- c(rep(1,0),seq(1,fs[i],length.out=niterny-0))
	yny[,i] <- as.integer(dPOSny*p)
	lines(doyny,yny[,i],lty=2,col=cols[i])
}
	mtext(side=2,line=2.5,'Daily Cases')
	mtext('New York')
plot(doywa,dPOSwa,ylim=c(0,1.5E3),type='l',lwd=2,bty='n',xlim=c(60,101))
fs <- seq(0.25,2.5,length.out=10)
ywa <- matrix(NA,ncol=length(fs),nrow=niterwa)
for(i in 1:length(fs)){
	p <- c(rep(1,0),seq(1,fs[i],length.out=niterwa-0))
	ywa[,i] <- as.integer(dPOSwa*p)
	lines(doywa,ywa[,i],lty=2,col=cols[i])
}
	mtext(side=2,line=2.5,'Daily Cases')
	mtext('Washington')

#####################################################################################
## ESTIMATE PARAMETERS ##############################################################
#####################################################################################
##--COMPILE STAN CODE--###############
mod_b   <- stan_model('stan_SIER_beta.stan') #compile the stan code
mod_b_s <- stan_model('stan_SIER_beta_sigma.stan') #compile the stan code

##--CODED WITH EULER FORWARD--####
mod     <- stan_model('stan_SEIR_discrete.stan')

##--TIME-DEPENDENT BETA--######
mod_betaAR1 <- stan_model('stan_SEIR_discrete_beta_AR1.stan')

#################################
## FIT TO NY DATA ###############
#################################
E0 <- 0.00001
S0 <- 1.0 - E0
R0 <- 0.0
I0 <- 0.0

x0 <- c(S0,E0,I0,R0) #initial conditions

MCMCny <- list()
for(i in 1:length(fs)){
print(i)
	d <- yny[,i] 
	data <- list(POP=POPny,
				 y=d,
				 N_obs=length(d),
				 t_obs=1:length(d),
				 x0 = x0)
	
	MCMCny[[i]] <- sampling(mod_b,data=data,open_progress=TRUE)
}

opt <- optimizing(mod_betaAR1,data=data,algorithm="Newton",hessian=TRUE)
plot(opt$par[1:38])
abline(h=1)
plot(opt$par[1:38])
lines(opt$par[1:38] + 2*sqrt(diag(solve(-opt$hessian)))[1:38])
lines(opt$par[1:38] - 2*sqrt(diag(solve(-opt$hessian)))[1:38])

opt <- optimizing(mod_betaAR1,data=data,algorithm="Newton",hessian=TRUE)

mcmc <- sampling(mod_betaAR1,data=data,open_progress=TRUE)
#################################
## FIT TO WA DATA ###############
#################################
E0 <- 0.00001
S0 <- 1.0 - E0
R0 <- 0.0
I0 <- 0.0

x0 <- c(S0,E0,I0,R0) #initial conditions

MCMCwa <- list()
for(i in 1:length(fs)){
print(i)
	d <- ywa[,i] 
	data <- list(POP=POPwa,
				 y=d,
				 N_obs=length(d),
				 t_obs=1:length(d),
				 x0 = x0)
	
	MCMCwa[[i]] <- sampling(mod_b,data=data,open_progress=TRUE)
}

####################################
## PLOTS ###########################
####################################
cols <- terrain.colors(15)

bw <- 0.01
par(mfrow=c(2,1))
post <- extract(MCMCny[[1]])
den <- density(post$theta/(1/7),bw=bw)
plot(den$x,den$y,type='l',ylim=c(0,60),xlim=c(2,4.2),col=cols[1],bty='n')
abline(v=mean(post$theta)/(1/7),col=cols[1],lty=2)
for(i in 2:length(fs)){
	post <- extract(MCMCny[[i]])
	den <- density(post$theta/(1/7),bw=bw)
	lines(den$x,den$y,col=cols[i])
	abline(v=mean(post$theta)/(1/7),col=cols[i],lty=2)
}
mtext('New York')

post <- extract(MCMCwa[[1]])
den <- density(post$theta/(1/7),bw=bw)
plot(den$x,den$y,type='l',ylim=c(0,60),xlim=c(2,4.2),col=cols[1],bty='n')
abline(v=mean(post$theta)/(1/7),col=cols[1],lty=2)
for(i in 2:length(fs)){
	post <- extract(MCMCwa[[i]])
	den <- density(post$theta/(1/7),bw=bw)
	lines(den$x,den$y,col=cols[i])
	abline(v=mean(post$theta)/(1/7),col=cols[1],lty=2)
}
mtext('Washington')
mtext(side=1,expression(italic('R'[0])),line=2.5,cex=1.25)

par(mfrow=c(1,1))
post <- extract(MCMCny[[1]])
hist(post$theta,xlim=c(0.48,0.60))
for(i in 2:length(fs)){
	post <- extract(MCMCny[[i]])
	hist(post$theta,add=TRUE,col=i)
}


par(mfrow=c(1,1))
post <- extract(MCMCwa[[1]])
hist(post$theta,xlim=c(0.48,0.60))
for(i in 2:length(fs)){
	post <- extract(MCMCwa[[i]])
	hist(post$theta,add=TRUE,col=i)
}













#sqrt(solve(-opt_SEIR$hessian))

opt_SEIR_new <- optimizing(mod_SEIR_new,data=data,hessian=TRUE,init='0')
opt_SEIR$par[1]/(1/14)


##--MCMC--######
mcmc_SEIR <- sampling(mod_SEIR,data=data,open_progress=TRUE) 
post_SEIR <- extract(mcmc_SEIR)

mcmc_SEIR_new <- sampling(mod_SEIR_new,data=data,open_progress=TRUE)

##--PLOTS--###
gamma <- 1/14
xin <- seq(0,15,0.01)

par(mfrow=c(1,2))
hist(post_SEIR$theta,main=''); mtext(expression(beta),cex=2)

plot(xin,dnorm(xin,mean=0.2/gamma,sd=3),type='l',bty='n')
par(new=TRUE)
hist(post_SEIR$theta/(1/14),freq=FALSE,breaks=10,xlim=c(0,15),xaxt='n',yaxt='n',main='')
	mtext(expression('Prior and Posterior R'[0]),cex=1.5)




