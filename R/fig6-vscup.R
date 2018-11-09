rm(list=ls(all=TRUE)) 
source("supercoev.R")
require(tikzDevice)
#############################################################################
############## alpha,gamma vs cup (up-regulation cost, figure 6a and 6b)

#### Draw directly the figure:
thedat <- read.table(file = "dat/fig6-vscup.csv", header=T)
# draw figure
epsoutput("figs/fig6-alpha.eps")
with(thedat,plot(cup,alphahi,type="l",lwd=2,ylim=c(0,4.5),xlab=expression("Up-regulation cost" ~ (c_up)),ylab=expression("CoESS virulence" ~ (alpha)), cex.lab=1.1))
with(subset(thedat,alpharep>0),lines(cup,alpharep,lwd=2,lty=2))
with(subset(thedat),lines(cup,alphalo,lwd=2,lty=1))
dev.off()

epsoutput("figs/fig6-gamma.eps")
with(thedat,plot(cup,gammahi,type="l",lwd=2,ylim=c(0,4.5),xlab=expression("Up-regulation cost" ~ (lambda)),ylab=expression("CoESS recovery rate" ~ (gamma)), cex.lab=1.1))
with(subset(thedat,alpharep>=0),lines(cup,gammarep,lwd=2,lty=2))
with(subset(thedat),lines(cup,gammalo,lwd=2,lty=1))
dev.off()

####################################
###### or run the simulations: #####

kappa<-0.01
betam<-10
mu<-1
rmax<-2
sigma0<-0.05
epsilon<-1
c<-0.05

thedat<-NULL
cup<-0.0
parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
# solve for cup=0
hicoess  <- multiroot(f = selgrad, start = c(3,4),positive=T) # upper coess
eqhi<-equilibrium(c((hicoess$root)[1],(hicoess$root)[2])) # equilibrium of susceptible and infected for the higher CoESS
repellor <- multiroot(f = selgrad, start = c(1,0),positive=T) # repellor 
locoess  <- uniroot(selgraduni, c(1,3), tol = 0.0001, maxiter = 1000) # loweer coess
eqlo<-equilibrium(c(locoess$root,0)) #  equilibrium of susceptible and infected for the lower CoESS

thedat<-rbind(thedat,data.frame(cup=0.0,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=(repellor$root)[1],gammarep=(repellor$root[2]),alphalo=locoess$root,gammalo=0,eqhiS=eqhi[1],eqhiI=eqhi[2],eqloS=eqlo[1],eqloI=eqlo[2]))

computeall=T
for (cup in seq(0.0,0.0346,by=0.0001)) {
  parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
  
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
  
  thedat<-rbind(thedat,data.frame(cup=cup,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=(repellor$root)[1],gammarep=(repellor$root[2]),alphalo=locoess$root,gammalo=0,eqhiS=NA,eqhiI=NA,eqloS=eqlo[1],eqloI=eqlo[2]))
}

for (cup in seq(0.039,0.05,by=0.001)) {
  
  parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
  
  locoess  <- uniroot(selgraduni, c(locoess$root-0.5,locoess$root+0.5), tol = 0.0001, maxiter = 1000)
  eqlo<-equilibrium(c(locoess$root,0))
  
  thedat<-rbind(thedat,data.frame(cup=cup,alphahi=NA,gammahi=NA,alpharep=NA,gammarep=NA,alphalo=locoess$root,gammalo=0,eqhiS=NA,eqhiI=NA,eqloS=eqlo[1],eqloI=eqlo[2]))
}
# save in csv table
write.table(thedat,"dat/fig6-vscup.csv")

# draw figure
epsoutput("figs/fig6-alpha.eps")
with(thedat,plot(cup,alphahi,type="l",lwd=2,ylim=c(0,4.5),xlab=expression("Up-regulation cost" ~ (c_up)),ylab=expression("CoESS virulence" ~ (alpha)), cex.lab=1.1))
with(subset(thedat,alpharep>0),lines(cup,alpharep,lwd=2,lty=2))
with(subset(thedat),lines(cup,alphalo,lwd=2,lty=1))
dev.off()

epsoutput("figs/fig6-gamma.eps")
with(thedat,plot(cup,gammahi,type="l",lwd=2,ylim=c(0,4.5),xlab=expression("Up-regulation cost" ~ (lambda)),ylab=expression("CoESS recovery rate" ~ (gamma)), cex.lab=1.1))
with(subset(thedat,alpharep>=0),lines(cup,gammarep,lwd=2,lty=2))
with(subset(thedat),lines(cup,gammalo,lwd=2,lty=1))
dev.off()

# only if needs the tex output
# tikz('figs/fig6-alpha.tex', width=3, height=3.5)
# with(thedat,plot(cup,alphahi,type="l",lwd=2,ylim=c(0,4.5),xlab="Up-regulation cost ($c_{up}$)",ylab="CoESS virulence ($\\alpha^*$)"))
# with(subset(thedat),lines(cup,alpharep,lwd=2,lty=2))
# with(subset(thedat),lines(cup,alphalo,lwd=2,lty=1))
# dev.off()      
# 
# tikz('figs/fig6-gamma.tex', width=3, height=3.5)
# with(thedat,plot(cup,gammahi,type="l",lwd=2,ylim=c(0,4.5),xlab="Up-regulation cost ($c_{up}$)",ylab="CoESS recovery rate ($\\gamma^*$)"))
# with(subset(thedat,alpharep>=0),lines(cup,gammarep,lwd=2,lty=2))
# with(subset(thedat),lines(cup,gammalo,lwd=2,lty=1))
# dev.off()     