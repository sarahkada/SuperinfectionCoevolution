rm(list=ls(all=TRUE))
source("supercoev.R")

#############################################################################
############## alpha,gamma vs sigma (figure 4)

### to draw directly the figure, open and plot:
thedat <- read.table(file = "dat/vssigma-fig3-kappa.csv", header=T)

epsoutput(paste("figs/fig3-alpha-kappa.eps"))
with(thedat,plot(sigma,alphahi,type="l",log="x",lwd=2,ylim=c(0,4.5),ylab="CoESS virulence"))
with(subset(thedat,alpharep>0),lines(sigma,alpharep,lwd=2,lty=2))
with(subset(thedat,alpharep>0),lines(sigma,alphalo,lwd=2,lty=1))
dev.off()

epsoutput(paste("figs/fig3-gamma-kappa.eps"))
with(thedat,plot(sigma,gammahi,type="l",log="x",lwd=2,ylim=c(0,4),ylab="CoESS virulence",xlab=expression(paste("Susceptibility to superinfection (",sigma,")",sep=""))))
with(subset(thedat,gammarep>=0),lines(sigma,gammarep,lwd=2,lty=2))
with(subset(thedat,gammarep>=0),lines(sigma,gammalo,lwd=2,lty=1))
dev.off()
####################################
# or run the simulations: #########
kappa<-0.01
betam<-10
mu<-1
rmax<-2
sigma0<-0
epsilon<-1
c<-0.05
cup<-0

thedat<-NULL
sigma0<-0.01
parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)

hicoess  <- multiroot(f = selgrad, start = c(3,4),positive=T)
eqhi<-equilibrium(c((hicoess$root)[1],(hicoess$root)[2]))
repellor <- multiroot(f = selgrad, start = c(1,0),positive=T)
locoess  <- uniroot(selgraduni, c(1,3), tol = 0.0001, maxiter = 1000)
eqlo<-equilibrium(c(locoess$root,0))

thedat<-rbind(thedat,data.frame(sigma= sigma0,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=(repellor$root)[1],gammarep=(repellor$root[2]),alphalo=locoess$root,gammalo=0,eqhiS=eqhi[1],eqhiI=eqhi[2],eqloS=eqlo[1],eqloI=eqlo[2]))

computeall=T

for (sig in c(seq(0.01,0.1,by=0.005),seq(0.1,1.0,by=0.01),seq(1,10))) {
  
  sigma0<-sig
  
  parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
  
  hicoess  <- multiroot(f = selgrad, start = c((hicoess$root)[1],(hicoess$root)[2]),positive=T)
  eqhi<-equilibrium(c((hicoess$root)[1],(hicoess$root)[2]))
  
  if(computeall) {
    repellor <- multiroot(f = selgrad, start = c((repellor$root)[1],(repellor$root)[2]),positive=T)
    locoess  <- uniroot(selgraduni, c(locoess$root-0.5,locoess$root+0.5), tol = 0.0001, maxiter = 1000)
    eqlo<-equilibrium(c(locoess$root,0))
  }
  else {
    repellor$root=c(-5,NA) 
    locoess$root=-1
    eqlo<-c(NA,NA)
  }
  
  # condition to stop computation of low CoESS and repellor when the repellor and low CoESS collide:
  if(abs((repellor$root)[1]-(locoess$root))<0.005) {
    computeall=F
    print('stop computation for low coess and repellor')
  }
  
  thedat<-rbind(thedat,data.frame(sigma=sig,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=(repellor$root)[1],gammarep=(repellor$root[2]),alphalo=locoess$root,gammalo=0,eqhiS=eqhi[1],eqhiI=eqhi[2],eqloS=eqlo[1],eqloI=eqlo[2]))
}

write.table(thedat,paste("dat/vssigma-fig3-kappa.csv"))

epsoutput(paste("figs/fig3-alpha-kappa.eps"))
with(thedat,plot(sigma,alphahi,type="l",log="x",lwd=2,ylim=c(0,4.5),ylab="CoESS virulence"))
with(subset(thedat,alpharep>0),lines(sigma,alpharep,lwd=2,lty=2))
with(subset(thedat,alpharep>0),lines(sigma,alphalo,lwd=2,lty=1))
dev.off()

epsoutput(paste("figs/fig3-gamma-kappa.eps"))
with(thedat,plot(sigma,gammahi,type="l",log="x",lwd=2,ylim=c(0,4),ylab="CoESS virulence",xlab=expression(paste("Susceptibility to superinfection (",sigma,")",sep=""))))
with(subset(thedat,gammarep>=0),lines(sigma,gammarep,lwd=2,lty=2))
with(subset(thedat,gammarep>=0),lines(sigma,gammalo,lwd=2,lty=1))
dev.off()

# # export files in .tex
# require(tikzDevice)
# tikz('figs/fig4-alpha-kappa.tex', width=3, height=3.5, standAlone=T)
# with(thedat,plot(sigma,alphahi,type="l",log="x",lwd=2,ylim=c(0.01,4.5),ylab="CoESS virulence ($\\alpha^*$)",xlab="Susceptibility to superinfection ($\\sigma$)"))
# with(subset(thedat,alpharep>0),lines(sigma,alpharep,lwd=2,lty=2))
# with(subset(thedat,alpharep>0),lines(sigma,alphalo,lwd=2,lty=1))
# dev.off()      

# tikz('figs/fig4-gamma-kappa.tex', width=3, height=3.5,standAlone=T)
# with(thedat,plot(sigma,gammahi,type="l",log="x",lwd=2,ylim=c(0,4),ylab="CoESS recovery rate ($\\gamma^*$)",xlab="Susceptibility to superinfection ($\\sigma$)"))
# with(subset(thedat,gammarep>=0),lines(sigma,gammarep,lwd=2,lty=2))
# with(subset(thedat,gammarep>=0),lines(sigma,gammalo,lwd=2,lty=1))
# dev.off()      