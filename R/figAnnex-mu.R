rm(list=ls(all=TRUE))
source("supercoev.R")
require(tikzDevice)
#############################################################################
############## alpha,gamma vs mu (Appendix figure A1-b)
# Attention !! Use very small step throughout the computation !!

#### Draw directly the figure:
thedat <- read.table(file = "dat/figAnnexe-mu.csv", header=T)

epsoutput("figs/figAnnexe-mu-alpha.eps")
with(thedat,plot(mu,alphahi,type="l",lwd=2,ylim=c(0,4),xlab=expression("Mortality" ~ (mu)),ylab=expression("CoESS virulence" ~ (alpha)), cex.lab=1.1))
with(subset(thedat,alpharep>0),lines(mu,alpharep,lwd=2,lty=2))
with(subset(thedat),lines(mu,alphalo,lwd=2,lty=1))
dev.off()

epsoutput("figs/figAnnexe-mu-gamma.eps")
with(thedat,plot(mu,gammahi,type="l",lwd=2,ylim=c(0,12),xlab=expression("Mortality" ~ (mu)),ylab=expression("CoESS recovery rate" ~ (gamma)), cex.lab=1.1))
with(subset(thedat,alpharep>=0),lines(mu,gammarep,lwd=2,lty=2))
with(subset(thedat),lines(mu,gammalo,lwd=2,lty=1))
dev.off()


############## or run the simulations: 

kappa<-0.01
betam<-10
c<-0.05
rmax<-2
epsilon<-1
cup<-0.0
sigma0<-0.05

thedat<-NULL 
mu<-0.5
parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
hicoess  <- multiroot(f = selgrad, start= c(4,10),positive=T)
repellor <- multiroot(f = selgrad, start = c(2,2.1),positive=T)
locoess  <- uniroot(selgraduni, c(1,2), tol = 0.0001, maxiter = 1000)
thedat<-rbind(thedat,data.frame(mu=mu,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=(repellor$root)[1],gammarep=(repellor$root[2]),alphalo=locoess$root,gammalo=0))

computeall=T

for (mu in c(seq(0.5, 1.07, by=0.0005))) {
  parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
  tryCatch(
{
  hicoess  <- multiroot(f = selgrad, start = c((hicoess$root)[1],(hicoess$root)[2]),positive=T)
  if(computeall) {
    repellor <- multiroot(f = selgrad, start = c((repellor$root)[1],(repellor$root)[2]),positive=T)
    locoess  <- uniroot(selgraduni, c(locoess$root-1.5,locoess$root+.5), tol = 0.0001, maxiter = 1000)
  }
  
  if(abs((repellor$root)[1]-(locoess$root))<0.01)
  {
    computeall=F
  }
  thedat<-rbind(thedat,data.frame(mu=mu,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=(repellor$root)[1],gammarep=(repellor$root[2]),alphalo=locoess$root,gammalo=0))
}
  )
}

mu<-1.07
parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
hicoess  <- multiroot(f = selgrad, start= c(2,2.9),positive=T)

for (mu in c(seq(1.07,1.2, by=0.0005))) {
  parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
  tryCatch(
{
  hicoess  <- multiroot(f = selgrad, start = c((hicoess$root)[1],(hicoess$root)[2]),positive=T)
  thedat<-rbind(thedat,data.frame(mu=mu,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=NA,gammarep=NA,alphalo=NA,gammalo=NA))
}
  )
}

write.table(thedat,paste("dat/figAnnexe-mu",".csv",sep=""))

epsoutput("figs/figAnnexe-mu-alpha.eps")
with(thedat,plot(mu,alphahi,type="l",lwd=2,ylim=c(0,4),xlab=expression("Mortality" ~ (mu)),ylab=expression("CoESS virulence" ~ (alpha)), cex.lab=1.1))
with(subset(thedat,alpharep>0),lines(mu,alpharep,lwd=2,lty=2))
with(subset(thedat),lines(mu,alphalo,lwd=2,lty=1))
dev.off()

epsoutput("figs/figAnnexe-mu-gamma.eps")
with(thedat,plot(mu,gammahi,type="l",lwd=2,ylim=c(0,12),xlab=expression("Mortality" ~ (mu)),ylab=expression("CoESS recovery rate" ~ (gamma)), cex.lab=1.1))
with(subset(thedat,alpharep>=0),lines(mu,gammarep,lwd=2,lty=2))
with(subset(thedat),lines(mu,gammalo,lwd=2,lty=1))
dev.off()

# # exporting the tex files:
# tikz('figs/figAnnexe-mu-alpha.tex', width=3, height=3.5, standAlone=TRUE)
# with(thedat,plot(mu,alphahi,type="l",lwd=2,ylim=c(0,4),xlab= "",ylab=""))
# with(subset(thedat,alpharep>0),lines(mu,alpharep,lwd=2,lty=2))
# with(subset(thedat),lines(mu,alphalo,lwd=2,lty=1))
# dev.off()      
# 
# tikz('figs/figAnnexe-mu-gamma.tex', width=3, height=3.5, standAlone=TRUE)
# with(thedat,plot(mu,gammahi,type="l",lwd=2,ylim=c(0,12),xlab="Natural mortality rate ($\\mu$)",ylab=""))
# with(subset(thedat,alpharep>=0),lines(mu,gammarep,lwd=2,lty=2))
# with(subset(thedat),lines(mu,gammalo,lwd=2,lty=1))
# dev.off()