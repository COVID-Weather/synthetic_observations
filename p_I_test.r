library(fields)


dt <- 0.1
niter <- 400/0.1
nsamp <- 200

Npop <- 1E4
I=S=E=R <- numeric(niter)
Iobs=Iobs_ptest <- numeric(niter)

pa <- c(0.001,0.002,0.005,0.01,0.01,0.05,0.1)
np <- length(pa)
ptest1 <- matrix(0,ncol=niter,nrow=np)
ptest2 <- matrix(0,ncol=niter,nrow=np)

E[1] <- 0.001
S[1] <- 1 - E[1]

gamma <- 0.2
sigma <- 0.1
beta <- 0.4

cols <- terrain.colors(np+3)

for(i in 2:niter){
	S[i] = S[i-1] - dt*(beta*I[i-1]*S[i-1])
	E[i] = E[i-1] + dt*(beta*I[i-1]*S[i-1] - sigma*E[i-1])
	I[i] = I[i-1] + dt*(sigma*E[i-1] - gamma*I[i-1])
	R[i] = R[i-1] + dt*(gamma*I[i-1])
	
	Iobs[i] <- rpois(1,lambda=I[i]*Npop)
	for(j in 1:np){
		ptest1[j,i] <-  I[i]/(I[i] + pa[j]*(S[i] + E[i] + R[i]))
		ptest2[j,i] <-  I[i]^2/(I[i]^2 + (pa[j]/500)*(S[i] + E[i] + R[i]))
	}
}


#####################################################################
## MULTIPLE BETAS ###################################################
#####################################################################
niter <- 400/dt
nbeta <- 5

Npop <- 1E4
I=S=E=R <- matrix(0,ncol=niter,nrow=nbeta)

Sobs=Eobs=Iobs=Robs <- matrix(0,ncol=niter,nrow=nsamp)

E[,1] <- 0.001
S[,1] <- 1 - E[,1]

gamma <- 1/7
sigma <- 1/5
beta  <- seq(1/5,0.5,length.out=nbeta)

for(i in 2:niter){
	S[,i] = S[,i-1] - dt*(beta*I[,i-1]*S[,i-1])
	E[,i] = E[,i-1] + dt*(beta*I[,i-1]*S[,i-1] - sigma*E[,i-1])
	I[,i] = I[,i-1] + dt*(sigma*E[,i-1] - gamma*I[,i-1])
	R[,i] = R[,i-1] + dt*(gamma*I[,i-1])
	
#	Sobs[,i] <- rpois(nsamp,lambda=S[,i]*Npop)
#	Eobs[,i] <- rpois(nsamp,lambda=E[,i]*Npop)
#	Iobs[,i] <- rpois(nsamp,lambda=I[,i]*Npop)
#	Robs[,i] <- rpois(nsamp,lambda=R[,i]*Npop)
}

par(mfrow=c(2,2))
matplot(t(S),type='l'); mtext('Susceptible')
matplot(t(E),type='l'); mtext('Exposed')
	legend('topright',legend=c(paste(round(beta/gamma,1))),bty='n',lty=1:5,col=1:5)
matplot(t(I),type='l'); mtext('Infected')
matplot(t(R),type='l'); mtext('Removed')

par(mfrow=c(1,2))
plot(-999,xlim=c(1,400),ylim=c(0,0.25))
	for(j in 1:nbeta){
		lines(I[j,]/(I[j,]+S[j,]+E[j,]+R[j,]),col=j,lty=j,lwd=1.35)
	}
mtext(side=2,'I/(I+S+E+R)',line=2.5)
plot(-999,xlim=c(1,400),ylim=c(0,0.25))
	for(j in 1:nbeta){
		lines(6*I[j,]^2/(6*I[j,]^2+S[j,]+E[j,]+R[j,]),col=j,lty=j,lwd=1.35)
	}
legend('topright',legend=c(paste(round(beta/gamma,1))),bty='n',lty=1:5,col=1:5)
mtext(side=2,expression('aI'^2*'/(aI'^2*'S+E+I+R)'),line=2.5)
mtext(side=1,'Days since Beginning of Infection',outer=TRUE)




logs <- seq(-5,0,length.out=5)
s    <- exp(logs)

par(mfrow=c(2,1))
plot(-999,xlim=c(1,400),ylim=c(0,1))
for(i in 1:length(s)){
	for(j in 1:nbeta){
		lines(I[j,]/(I[j,]+s[i]*(S[j,]+E[j,]+R[j,])),col=j,lty=i)
	}
}
mtext(side=2,expression('I/(I+s(S+E+R))'),line=2.5)
plot(-999,xlim=c(1,400),ylim=c(0,1))
for(i in 1:length(s)){
	for(j in 1:nbeta){
		lines(I[j,]^2/(I[j,]^2+(s[i]/10)*(S[j,]+E[j,]+R[j,])),col=j,lty=i)
	}
}
mtext(side=2,expression('I'^2*'/(I'^2*'+s(S+E+R))'),line=2.5)
legend(x=290,y=1.05,bty='n',legend=paste(round(s,3)),lty=1:length(s),col=1)
legend(x=190,y=1.05,bty='n',legend=paste(round(beta/gamma,1)),lty=1,col=1:length(s))
mtext(outer=TRUE,side=1,'Days since Beginning of Infection')






########################################################
# par(mfrow=c(1,2))
# plot(-999,xlim=c(0,niter),ylim=c(0,1),bty='n')
# for(i in 1:np){
	# lines(ptest1[i,],col=cols[i])
# }
# legend('topright',legend=paste(round(pa,5)),cex=0.4,lty=1,col=cols[1:np],bty='n')
# plot(-999,xlim=c(0,niter),ylim=c(0,1),bty='n')
# for(i in 1:np){
	# lines(ptest2[i,],col=cols[i])
# }


#########################################################

sts <- c('NY','FL','WA','TX','VA','IL')
pdf('d:/dropbox/working/covid19/plots/p_I_test_I_SEIR_trajectories_from_prior_STATES.pdf',height=7,width=9)
par(mfrow=c(3,2),mar=c(2,2,2,2),oma=c(2,2,2,2))
for(p in 1:length(sts)){
	plot(DATA[[sts[p]]]$dPOS,ylim=c(0,max(DATA[[sts[p]]]$Ntest)*0.4),pch=19,type='b',bty='n',xlim=c(0,40))
	for(i in 1:length(s)){
		for(j in 1:nbeta){
			ntest <- DATA[[sts[p]]]$Ntest
			ptest <- I[j,]/(I[j,]+s[i]*(S[j,]+E[j,]+R[j,]))
			lines(ntest*ptest[1:length(ntest)],col=adjustcolor(j,alpha.f=0.6),lty=i)
		}
	}
	mtext(sts[p],adj=0.15,line=-2)
}
mtext(side=1,outer=TRUE,'Days since Beginning of Infection',line=0.5)
mtext(side=2,outer=TRUE,'Daily New Infections')
dev.off()

##########################################################

logss <- seq(-10,0,length.out=10)
ss    <- exp(logss)

xi <- 100/dt
xii <- seq(0,100,length.out=1000)

cols <- terrain.colors(13)

par(mfrow=c(2,2),mar=c(2,2,2,4),oma=c(2,2,2,2))
for(p in c(1,2,3,5)){
plot(-999,xlim=c(0,100),ylim=c(0,1))
	for(i in 1:length(ss)){
		lines(xii,I[p,1:xi]/(I[p,1:xi]+ss[i]*(S[p,1:xi]+E[p,1:xi]+R[p,1:xi])),col=cols[i])
	}
	lines(xii,I[p,1:xi]/(I[p,1:xi]+ss[i]*(S[p,1:xi]+E[p,1:xi]+R[p,1:xi])),col='black',lwd=2)
	mtext(beta[p]/gamma)
}
image.plot(matrix(ss),legend.only=TRUE,col=terrain.colors(100)[1:80],outer=TRUE)
mtext(side=1,outer=TRUE,'Days since Beginning of Infection',line=0.4)
mtext(side=2,outer=TRUE,'I/(I+s(S+E+R))',line=0.4)





#########################################################
fs <- seq(0.25,2.5,length.out=np)
cols <- terrain.colors(12)

Nout <- 1:200

par(mfrow=c(3,1),mar=c(2,2,2,2),oma=c(2,2,2,2))
plot(-999,xlim=range(Nout),ylim=c(0,3),bty='n')
for(i in 1:length(fs)){
	p <- c(rep(1,0),seq(1,fs[i],length.out=length(Nout)-0))
	lines(Nout,p,lty=2,col=cols[i])
	abline(h=1)
}
mtext(side=2,line=2.5,'Case Multiplier')
plot(-999,xlim=range(Nout),ylim=c(0,3),bty='n')
for(i in 1:length(fs)){
	p <- c(rep(1,0),seq(1,fs[i],length.out=length(Nout)-0))
	lines(Nout,p*ptest[i,1:200]+p,lty=2,col=cols[i])
	abline(h=1)
}




