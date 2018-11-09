rm(list=ls(all=TRUE))
source("supercoev.R")
require(tikzDevice)
#############################################################################
############## alpha,gamma vs c (Appendix)

#### Draw directly the figure:
thedat <- read.table(file = "dat/figannexe-c.csv", header=T)

epsoutput("figs/figAnnex-c-alpha.eps")
with(thedat,plot(c,alphahi,type="l",lwd=2,ylim=c(0,6),xlab=expression("Up-regulation cost" ~ (c_up)),ylab=expression("CoESS virulence" ~ (alpha)), cex.lab=1.1))
with(subset(thedat,alpharep>0),lines(c,alpharep,lwd=2,lty=2))
with(subset(thedat),lines(c,alphalo,lwd=2,lty=1))
dev.off()

epsoutput("figs/figAnnex-c-gamma.eps")
with(thedat,plot(c,gammahi,type="l",lwd=2,ylim=c(0,33),xlab=expression("Cost of immune system maintenantce" ~ (c)),ylab=expression("CoESS recovery rate" ~ (gamma)), cex.lab=1.1))
with(subset(thedat,alpharep>=0),lines(c,gammarep,lwd=2,lty=2))
with(subset(thedat),lines(c,gammalo,lwd=2,lty=1))
dev.off()

####################################
# or run the simulations: #########
kappa<-0.01
betam<-10
mu<-1
sigma0<-0.05
epsilon<-1
cup<-0.0
rmax<-2

thedat<- NULL

c<-0.015
parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
hicoess  <- multiroot(f = selgrad, start = c(4.5,19.7),positive=T)
thedat<-rbind(thedat,data.frame(c=c,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=NA,gammarep=NA,alphalo=NA,gammalo=NA))
computeall=T
for (c in seq(0.015,0.019,by=0.0005)) {  
  tryCatch(
{
  parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
  
  if(computeall) {
    hicoess  <- multiroot(f = selgrad, start = c((hicoess$root)[1],(hicoess$root)[2]),positive=T)
    # eqhi<-equilibrium(c((hicoess$root)[1],(hicoess$root)[2]))
  } else {
    hicoess$root=c(-1,NA)
  }
  thedat<-rbind(thedat,data.frame(c=c,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=NA,gammarep=NA,alphalo=NA,gammalo=NA))
}
  )
}

c<-0.019
parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
hicoess  <- multiroot(f = selgrad, start = c(4,14),positive=T)
thedat<-rbind(thedat,data.frame(c=c,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=NA,gammarep=NA,alphalo=NA,gammalo=NA))
for (c in seq(0.0200,0.046,by=0.001)) {  
  tryCatch(
{
  parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
  
  if(computeall) {
    hicoess  <- multiroot(f = selgrad, start = c((hicoess$root)[1],(hicoess$root)[2]),positive=T)
    # eqhi<-equilibrium(c((hicoess$root)[1],(hicoess$root)[2]))
  } else {
    hicoess$root=c(-1,NA)
  }
  thedat<-rbind(thedat,data.frame(c=c,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=NA,gammarep=NA,alphalo=NA,gammalo=NA))
}
  )
}

c<-0.046
parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
hicoess  <- multiroot(f = selgrad, start = c(2,4),positive=T)
repellor <- multiroot(f = selgrad, start = c(1.1,0),positive=T)
locoess  <- uniroot(selgraduni, c(0.0,3), tol = 0.0001, maxiter = 1000)
#eqlo<-equilibrium(c(locoess$root,0))
thedat<-rbind(thedat,data.frame(c=c,alphahi=NA,gammahi=NA,alpharep=(repellor$root)[1],gammarep=(repellor$root[2]),alphalo=locoess$root,gammalo=0))

#for (c in seq(0.046,0.061,by=0.0005)) { 
for (c in seq(0.046,0.0595,by=0.0005)) { 
  tryCatch(
{
  parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
  
  locoess  <- uniroot(selgraduni, c(locoess$root-0.5,locoess$root+0.5), tol = 0.0001, maxiter = 1000)
  print(c)
  if(computeall) {
    repellor <- multiroot(f = selgrad, start = c((repellor$root)[1],(repellor$root)[2]),positive=T)
    hicoess  <- multiroot(f = selgrad, start = c((hicoess$root)[1],(hicoess$root)[2]),positive=T)
  } else {
    repellor$root=c(-5,NA)
    hicoess$root=c(-1,NA)
  }
  if(abs((repellor$root)[1]-(hicoess$root)[1])<0.01)
  {
    computeall=F
  }
  
  thedat<-rbind(thedat,data.frame(c=c, alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=(repellor$root)[1],gammarep=(repellor$root[2]),alphalo=locoess$root,gammalo=0))
}
  )
}
#sols inf seule pour 0.06 et +
for (c in seq(0.061, 0.07,by=0.0005)) { 
  tryCatch(
{
  parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
  locoess  <- uniroot(selgraduni, c(locoess$root-0.5,locoess$root+0.5), tol = 0.0001, maxiter = 1000)
  thedat<-rbind(thedat,data.frame(c=c, alphahi=NA,gammahi=NA,alpharep=NA,gammarep=NA,alphalo=locoess$root,gammalo=0))
}
  )
}

# save thedat
write.table(thedat,"dat/figannexe-c.csv")

# draw the CoESS figure
epsoutput("figs/figAnnex-c-alpha.eps")
with(thedat,plot(c,alphahi,type="l",lwd=2,ylim=c(0,6),xlab=expression("Up-regulation cost" ~ (c_up)),ylab=expression("CoESS virulence" ~ (alpha)), cex.lab=1.1))
with(subset(thedat,alpharep>0),lines(c,alpharep,lwd=2,lty=2))
with(subset(thedat),lines(c,alphalo,lwd=2,lty=1))
dev.off()

epsoutput("figs/figAnnex-c-gamma.eps")
with(thedat,plot(c,gammahi,type="l",lwd=2,ylim=c(0,33),xlab=expression("Cost of immune system maintenantce" ~ (c)),ylab=expression("CoESS recovery rate" ~ (gamma)), cex.lab=1.1))
with(subset(thedat,alpharep>=0),lines(c,gammarep,lwd=2,lty=2))
with(subset(thedat),lines(c,gammalo,lwd=2,lty=1))
dev.off()

# if figure are needed in tex format:

# tikz('figs/figAnnex-c-alpha.tex', width=3, height=3.5)
# with(thedat,plot(c,alphahi,type="l",lwd=2,ylim=c(0,4), xlab="", ylab=""))
# with(subset(thedat,alpharep>0),lines(c,alpharep,lwd=2,lty=2))
# with(subset(thedat),lines(c,alphalo,lwd=2,lty=1))
# dev.off()

# tikz('figs/figAnnex-c-gamma.tex', width=3, height=3.5)
# with(thedat,plot(c,gammahi,type="l",lwd=2,ylim=c(0,19.5),xlab="Cost of defence maintenance ($c$)", ylab=""))
# with(subset(thedat,alpharep>=0),lines(c,gammarep,lwd=2,lty=2))
# with(subset(thedat),lines(c,gammalo,lwd=2,lty=1))
# dev.off()      