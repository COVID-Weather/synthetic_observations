
dat <- read.csv('d:/dropbox/working/covid19/states_daily.csv')

states <- unique(dat$state)

i=1

NY <- dat[dat$state==states[i],]

day   <- rev(NY$date)
dTEST <- rev(NY$totalTestResultsIncrease)
dPOS  <- rev(NY$positiveIncrease)
dNEG  <- rev(NY$negativeIncrease)

par(mfrow=c(2,2))
plot(day,dTEST); mtext(states[i])
plot(day,dPOS)
plot(day,dNEG)
plot(day,dPOS/dTEST,ylim=c(0,1))

i=i+1


plot(NY$date,NY$positive)

plot(NY$date,NY$hospitalized)


I0 <- 0.001
S0 <- 1.0 - I0
R0 <- 0.0
E0 <- 0.0
x0 <- c(S0,E0,I0,R0) #initial conditions

data <- list(POP=1E6,
		     y=cbind(NY$positive,NY$positive,NY$positive,NY$positive),
			 N_obs=length(NY$positive),
			 t_obs=1:length(NY$positive),
			 x0 = x0,
			 x_p=4,
			 N_obsvar=1,
			 theta_p=3,
			 i_obsvar=c(3))
opt_SEIR <- optimizing(mod_SEIR,data=data)
opt_SEIR$par[1:3]
opt_SEIR$par[2]/opt_SEIR$par[1]
