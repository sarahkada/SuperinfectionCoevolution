rm(list=ls(all=TRUE))
source("supercoev.R")

# debug tool: use traceback() to get where the error comes from
# use browser() inside a function to access its calculation
#############################################################################
############## alpha,gamma vs epsilon (fig 5a and 5b)

### to draw directly the figure, open and plot:
thedat <- read.table(file = "dat/fig5-vsepsilon.csv", header=T)

epsoutput("figs/fig5-alpha.eps")
with(thedat,plot(epsilon,alphahi,type="l",lwd=2,ylim=c(0,4.5),xlab=expression("Relative infected hosts fecundity"~ (epsilon)), ylab=expression("CoESS virulence" ~ (alpha)), cex.lab=1.1))
with(subset(thedat,alpharep>0),lines(epsilon,alpharep,lwd=2,lty=2))
with(subset(thedat),lines(epsilon,alphalo,lwd=2,lty=1))
dev.off()

epsoutput("figs/fig5-gamma.eps")
with(thedat,plot(epsilon,gammahi,type="l",lwd=2,ylim=c(0,6),xlab=expression("Relative infected hosts fecundity"~ (epsilon)), ylab=expression("CoESS recovery rate" ~ (gamma)), cex.lab=1.1))
with(subset(thedat,alpharep>=0),lines(epsilon,gammarep,lwd=2,lty=2))
with(subset(thedat),lines(epsilon,gammalo,lwd=2,lty=1))
dev.off()

####################################
# or run the simulations: #########
kappa<-0.01
betam<-10
mu<-1
rmax<-2
sigma0<-0.05
cup<-0.0
c<-0.05
thedat<-NULL

epsilon<-0.0
parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
hicoess  <- multiroot(f = selgrad, start = c(2,5),positive=T)
thedat<-rbind(thedat,data.frame(epsilon=epsilon,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=NA,gammarep=NA,alphalo=NA,gammalo=NA))

computeall=T

# Solve for higher coess only (from 0 to 0.96)
for (epsilon in seq(0, 0.97, by=0.01)) {
  if(computeall) {
    hicoess  <- multiroot(f = selgrad, start = c((hicoess$root)[1],(hicoess$root)[2]),positive=T)
  } else {
    hicoess$root=c(-1,NA)
  }
  thedat<-rbind(thedat,data.frame(epsilon=epsilon,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=NA,gammarep=NA,alphalo=NA,gammalo=NA))
}

epsilon<-0.97
parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
repellor <- multiroot(f = selgrad, start = c(1,0),positive=T)
locoess  <- uniroot(selgraduni, c(0.8,1.2), tol = 0.0001, maxiter = 1000)
thedat<-rbind(thedat,data.frame(epsilon=epsilon,alphahi=NA,gammahi=NA,alpharep=(repellor$root)[1],gammarep=(repellor$root[2]),alphalo=locoess$root,gammalo=0))

# repellor and hicoess sol from 0.97 to the end
  for (epsilon in seq(0.97, 1.093, by=0.005)) {
  locoess  <- uniroot(selgraduni, c(locoess$root-0.5,locoess$root+0.5), tol = 0.0001, maxiter = 1000)
  
  if(computeall) {
    repellor <- multiroot(f = selgrad, start = c((repellor$root)[1],(repellor$root)[2]),positive=T)
    hicoess  <- multiroot(f = selgrad, start = c((hicoess$root)[1],(hicoess$root)[2]),positive=T)
  } else {
    repellor$root=c(-5,NA)
    hicoess$root=c(-1,NA)
  }
  if(abs((repellor$root)[1]-(hicoess$root)[1])<0.005)
  {
    computeall=F
  }
  # thedat<-rbind(thedat,data.frame(epsilon=0.0,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=NA,gammarep=NA,alphalo=NA,gammalo=NA))
  thedat<-rbind(thedat,data.frame(epsilon=epsilon,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=(repellor$root)[1],gammarep=(repellor$root[2]),alphalo=locoess$root,gammalo=0))
}

# les solutions à partir de 1.10
for (epsilon in seq(1.10, 1.2, by=0.01)) {
  locoess  <- uniroot(selgraduni, c(locoess$root-0.5,locoess$root+0.5), tol = 0.0001, maxiter = 1000)
  # thedat<-rbind(thedat,data.frame(epsilon=0.0,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=NA,gammarep=NA,alphalo=NA,gammalo=NA))
  thedat<-rbind(thedat,data.frame(epsilon=epsilon,alphahi=NA,gammahi=NA,alpharep=NA,gammarep=NA,alphalo=locoess$root,gammalo=0))
}

write.table(thedat,"dat/fig5-vsepsilon.csv") # save new simulations
# plot:
epsoutput("figs/fig5-alpha.eps")
with(thedat,plot(epsilon,alphahi,type="l",lwd=2,ylim=c(0,4.5),xlab=expression("Relative infected hosts fecundity"~ (epsilon)), ylab=expression("CoESS virulence" ~ (alpha)), cex.lab=1.1))
with(subset(thedat,alpharep>0),lines(epsilon,alpharep,lwd=2,lty=2))
with(subset(thedat),lines(epsilon,alphalo,lwd=2,lty=1))
dev.off()

epsoutput("figs/fig5-gamma.eps")
with(thedat,plot(epsilon,gammahi,type="l",lwd=2,ylim=c(0,6),xlab=expression("Relative infected hosts fecundity"~ (epsilon)), ylab=expression("CoESS recovery rate" ~ (gamma)), cex.lab=1.1))
with(subset(thedat,alpharep>=0),lines(epsilon,gammarep,lwd=2,lty=2))
with(subset(thedat),lines(epsilon,gammalo,lwd=2,lty=1))
dev.off()

# if figure are needed in tex format:
# require(tikzDevice)
# tikz('figs/fig5-alpha.tex', width=3, height=3.5)
# with(thedat,plot(epsilon,alphahi,type="l",lwd=2,ylim=c(0,4.5),xlab="Relative infected hosts fecundity ($\\varepsilon$)", ylab="CoESS virulence ($\\alpha^*$)"))
# with(subset(thedat,alpharep>0),lines(epsilon,alpharep,lwd=2,lty=2))
# with(subset(thedat),lines(epsilon,alphalo,lwd=2,lty=1))
# dev.off()      
# 
# tikz('figs/fig5-gamma.tex', width=3, height=3.5)
# with(thedat,plot(epsilon,gammahi,type="l",lwd=2,ylim=c(0,6),xlab="Relative infected hosts fecundity ($\\varepsilon$)", ylab="CoESS recovery rate ($\\gamma^*$)"))
# with(subset(thedat,alpharep>=0),lines(epsilon,gammarep,lwd=2,lty=2))
# with(subset(thedat),lines(epsilon,gammalo,lwd=2,lty=1))
# dev.off()     