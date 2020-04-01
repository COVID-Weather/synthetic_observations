mod_SEIR <- stan_model('stan_SIER.stan') #compile the stan code

mod_SEIR_discrete <- stan_model('stan_SIER_discrete.stan') #compile the stan code

dat <- read.csv('d:/dropbox/working/covid19/south_korea.csv',stringsAsFactors=FALSE)
dat$doy <- as.numeric(strftime(as.character(paste(dat$year,dat$month,dat$day,sep='-')), format = "%j")) + dat$decimal_day

dat2 <- dat[10:nrow(dat),]

par(mfrow=c(2,2))
plot(dat$doy,dat$test)
plot(dat$doy,dat$positive)
plot(dat$doy,dat$suspected)
plot(dat$doy,dat$positive-dat$discharged)

bias <- 1/(dat$tests/dat$suspected)
bias[is.na(bias)] <- 1.0

plot(dat$positive/dat$test,ylim=c(0,0.1))


par(mfrow=c(3,2),mar=c(2,2,2,2))
plot(dat$doy,1/bias)
plot(dat$doy[-1],diff(dat$positive),ylim=c(0,1000),type='l')
lines(dat$doy[-1],bias[-1]*diff(dat$positive),ylim=c(0,1000),col='red')

plot(dat$doy,dat$positive,ylim=c(0,15000))
lines(dat$doy[-1],cumsum(bias[-1]*diff(dat$positive)))


plot(dat$doy,dat$positive - dat$discharged,ylim=c(0,15000))
lines(dat$doy[-1],cumsum(bias[-1]*diff(dat$positive)) - dat$discharged[-1])

plot(dat$doy,dat$discharged + dat$death)

N <- nrow(dat)
N <- nrow(dat2)

bias2 <- 1/(dat2$tests/dat2$suspected)
bias2[is.na(bias2)] <- 1.0


data_correct <- list(POP=2E4,
		     y=cbind(rep(1,N-1),rep(1,N-1), as.integer(cumsum(bias2[-1]*diff(dat2$positive)) - dat2$discharged[-1]) ,dat2$discharged[-1]),
			 N_obs=N-1,
			 t_obs=dat2$doy[-1],
			 x0 = x0,
			 x_p=4,
			 N_obsvar=2,
			 theta_p=3,
			 i_obsvar=c(3,4))

data <- list(POP=2E4,
		     y=cbind(rep(1,N),rep(1,N), dat2$positive - dat2$discharged, dat2$discharged),
			 N_obs=N,
			 t_obs=dat2$doy,
			 x0 = x0,
			 x_p=4,
			 N_obsvar=2,
			 theta_p=3,
			 i_obsvar=c(3,4))

opt_SEIR <- optimizing(mod_SEIR,data=data,init='0')
opt_SEIR$par[1:4]
opt_SEIR$par[1]/0.1

opt_SEIR <- optimizing(mod_SEIR,data=data_correct,init='0')
opt_SEIR$par[1:3]
opt_SEIR$par[1]/0.1


mcmc_SEIR <- sampling(mod_SEIR,data=data,open_progress=TRUE,init='0')


mcmc_SEIR_discrete <- sampling(mod_SEIR_discrete,data=data,open_progress=TRUE,init='0')


