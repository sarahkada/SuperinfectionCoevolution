rm(list=ls(all=TRUE))
source("supercoev-facilitation.R")
#############################################################################
############## alpha,gamma vs a
require(tikzDevice)

# only draw the figure:
thedat <- read.table(file = "dat/fig8-facilitation-table.csv", header=T)

epsoutput("figs/fig8-facilitation.eps")
with(thedat,plot(a,alphahi,type="l",lwd=2,xlim=c(-1.6,3), ylim=c(0,8),xlab=expression("a"),ylab=expression(alpha * "  " * gamma * "  "  * sigma), cex.lab=1.1))
with(subset(thedat),lines(a,gammahi,lwd=2,lty=2))
with(subset(thedat),lines(a,sigmavar,lwd=2,lty=3))
dev.off()

# re-run the simulations:
kappa<-0.01
betam<-10
mu<-1
rmax<- 2
epsilon<-1
c<-0.05
cup<-0.0
sigmaz<- 1.0
x<- 1
thedat<-NULL

a<- -1.6
parms<-c(a=a, kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigmaz=sigmaz, x=x, epsilon=epsilon, c=c, cup=cup)
hicoess  <- multiroot(f = selgrad, start = c(2,3),positive=T)
gammahi=(hicoess$root)[2]
sigmavar<- sigma(gammahi)
thedat<-rbind(thedat,data.frame(a=a,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2], sigmavar=sigmavar))

for (a in seq(-1.4,-1.01,by=0.01)) {
  parms<-c(a=a, kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigmaz=sigmaz, x=x, epsilon=epsilon, c=c, cup=cup)
  hicoess  <- multiroot(f = selgrad, start = c((hicoess$root)[1],(hicoess$root)[2]),positive=T)
  gammahi=(hicoess$root)[2]
  sigmavar<- sigma(gammahi)
  thedat<-rbind(thedat,data.frame(a=a,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2], sigmavar=sigmavar))
}

a<- -1.0
parms<-c(a=a, kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigmaz=sigmaz, x=x, epsilon=epsilon, c=c, cup=cup)
hicoess  <- multiroot(f = selgrad, start = c(2,3),positive=T)
gammahi=(hicoess$root)[2]
sigmavar<- sigma(gammahi)
thedat<-rbind(thedat,data.frame(a=a,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2], sigmavar=sigmavar))

for (a in seq(-1.0,-0.01,by=0.01)) {
  
  parms<-c(a=a, kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigmaz=sigmaz, x=x, epsilon=epsilon, c=c, cup=cup)
  
  hicoess  <- multiroot(f = selgrad, start = c((hicoess$root)[1],(hicoess$root)[2]),positive=T)
  gammahi=(hicoess$root)[2]
  sigmavar<- sigma(gammahi)
  thedat<-rbind(thedat,data.frame(a=a,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2], sigmavar=sigmavar))
}


a<- 0.0
parms<-c(a=a, kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigmaz=sigmaz, x=x, epsilon=epsilon, c=c, cup=cup)
hicoess  <- multiroot(f = selgrad, start = c(2,3),positive=T)
gammahi=(hicoess$root)[2]
sigmavar<- sigma(gammahi)
thedat<-rbind(thedat,data.frame(a=a,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2], sigmavar=sigmavar))

computeall=T
for (a in seq(0,3,by=0.01)) {
  
  parms<-c(a=a, kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigmaz=sigmaz, x=x, epsilon=epsilon, c=c, cup=cup)
  
  hicoess  <- multiroot(f = selgrad, start = c((hicoess$root)[1],(hicoess$root)[2]),positive=T)
  gammahi=(hicoess$root)[2]
  sigmavar<- sigma(gammahi)
  thedat<-rbind(thedat,data.frame(a=a,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2], sigmavar=sigmavar))
}

write.table(thedat,"dat/fig8-facilitation-table.csv")

epsoutput("figs/fig8-facilitation.eps")
with(thedat,plot(a,alphahi,type="l",lwd=2,xlim=c(-1.6,3), ylim=c(0,8),xlab=expression("a"),ylab=expression(alpha * "  " * gamma * "  "  * sigma), cex.lab=1.1))
with(subset(thedat),lines(a,gammahi,lwd=2,lty=2))
with(subset(thedat),lines(a,sigmavar,lwd=2,lty=3))
dev.off()

# only if needs the tex output
# tikz('figs/fig8-facilitation.tex', width=3, height=3.5, standAlone=T)
# with(thedat,plot(a,alphahi,type="l",lwd=2,xlim=c(-1.6,2.0), ylim=c(0,8),ylab="$\\alpha^* \\quad \\gamma^* \\quad \\sigma$",xlab="a"))
# with(subset(thedat),lines(a,gammahi,lwd=2,lty=2))
# with(subset(thedat),lines(a,sigmavar,lwd=2,lty=3))
# legend(-1.65, 8, c("$\\alpha^*$", "$\\gamma^*$", "$\\sigma$" ),lty=c(1, 2, 3), bty='y', cex=.75)
# dev.off()            
# 
# tikz('figs/fig8-facilitation1.tex', width=3, height=3.5, standAlone=T)
# with(thedat,plot(a,alphahi,type="l",lwd=2,xlim=c(-1.6,2.0), ylim=c(0,8),ylab="$\\alpha^* \\quad \\gamma^*$",xlab="a"))
# with(subset(thedat),lines(a,gammahi,lwd=2,lty=2))
# #with(subset(thedat),lines(a,sigmavar,lwd=2,lty=3))
# legend(-1.65, 8, c("$\\alpha^*$", "$\\gamma^*$" ),lty=c(1, 2), bty='y', cex=.75)
# dev.off()            
# 
# tikz('figs/fig8-facilitation2.tex', width=3, height=3.5, standAlone=T)
# with(thedat,plot(a,sigmavar,type="l",lwd=2,lty=3, xlim=c(-1.6,3.0), ylim=c(0,30),ylab="$\\sigma (\\gamma)$",xlab="a"))
# with(subset(thedat),lines(a,sigmavar,lwd=2,lty=3))
# #legend(-1.65, 8, c("$\\alpha^*$", "$\\gamma^*$", "$\\sigma$" ),lty=c( 3), bty='y', cex=.75)
# dev.off()            


