rm(list=ls(all=TRUE))
source("supercoev.R")
require(tikzDevice)
#############################################################################
############## alpha,gamma vs rm (Appendix figure A1-a)

#### Draw directly the figure:
thedat <- read.table(file = "dat/figAnnexe-rm.csv", header=T)

epsoutput("figs/figAnnex-rm-alpha.eps")
with(thedat,plot(rmax,alphahi,type="l",lwd=2,ylim=c(0,3.5),xlab=expression("Maximal reproduction rate" ~ (r_m)),ylab=expression("CoESS virulence" ~ (alpha)), cex.lab=1.1))
with(subset(thedat,alpharep>0),lines(rmax,alpharep,lwd=2,lty=2))
with(subset(thedat),lines(rmax,alphalo,lwd=2,lty=1))
dev.off()

epsoutput("figs/figAnnex-rm-gamma.eps")
with(thedat,plot(rmax,gammahi,type="l",lwd=2,ylim=c(0,4),xlim=c(1.2,2.4),xlab=expression("Maximal reproduction rate" ~ (r_m)),ylab=expression("CoESS recovery rate" ~ (gamma)), cex.lab=1.1))
with(subset(thedat,alpharep>=0),lines(rmax,gammarep,lwd=2,lty=2))
with(subset(thedat),lines(rmax,gammalo,lwd=2,lty=1))
dev.off()

####################################
# or run the simulations: #########
kappa<-0.01
betam<-10
mu<-1
sigma0<-0.05
epsilon<-1
cup<-0.0
c<-0.05

# sols pour COESS haute de 1.2 à 1.9.
thedat<-NULL
rmax<-1.2
parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
hicoess  <- multiroot(f = selgrad, start = c(1,0),positive=T)
thedat<-rbind(thedat,data.frame(rmax=rmax,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=NA,gammarep=NA,alphalo=NA,gammalo=NA))

computeall=T
for (rmax in seq(1.2,1.93,by=0.01)) {
  tryCatch(
{
  parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
  
  if(computeall) {
    hicoess  <- multiroot(f = selgrad, start = c((hicoess$root)[1],(hicoess$root)[2]),positive=T)
  } else {
    hicoess$root=c(1,NA)
  }
  
  thedat<-rbind(thedat,data.frame(rmax=rmax,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=NA,gammarep=NA,alphalo=NA,gammalo=NA))
},
error=function(cond) {
  message(cond)
  message(cat(" | Failed with g_step = ", g_step, "\n", sep=""))#, file=f))
},
warning=function(cond) {
  message(cond)
},
finally={
}
  )    
}

#sols repellor et inf for 1.9 to 2.3
rmax<-1.93
parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
repellor <- multiroot(f = selgrad, start = c(1.14,0),positive=T)
locoess  <- uniroot(selgraduni, c(0.0,5), tol = 0.0001, maxiter = 1000)
thedat <- rbind(thedat,data.frame(rmax=rmax,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=(repellor$root)[1],gammarep=(repellor$root[2]),alphalo=locoess$root,gammalo=0))

for (rmax in seq(1.93,2.29,by=0.005)) {
  
  parms<-c(kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
  
  locoess  <- uniroot(selgraduni, c(locoess$root-0.5,locoess$root+0.5), tol = 0.0001, maxiter = 1000)
  print(rmax)
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
  
  thedat<-rbind(thedat,data.frame(rmax=rmax, alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=(repellor$root)[1],gammarep=(repellor$root[2]),alphalo=locoess$root,gammalo=0))
}

write.table(thedat,"dat/figAnnexe-rm.csv")  # save table

epsoutput("figs/figAnnex-rm-alpha.eps") # plot CoESS virulence for varying rm
with(thedat,plot(rmax,alphahi,type="l",lwd=2,ylim=c(0,3.5),xlab=expression("Maximal reproduction rate" ~ (r_m)),ylab=expression("CoESS virulence" ~ (alpha)), cex.lab=1.1))
with(subset(thedat,alpharep>0),lines(rmax,alpharep,lwd=2,lty=2))
with(subset(thedat),lines(rmax,alphalo,lwd=2,lty=1))
dev.off()

epsoutput("figs/figAnnex-rm-gamma.eps")
with(thedat,plot(rmax,gammahi,type="l",lwd=2,ylim=c(0,4),xlim=c(1.2,2.4),xlab=expression("Maximal reproduction rate" ~ (r_m)),ylab=expression("CoESS recovery rate" ~ (gamma)), cex.lab=1.1))
with(subset(thedat,alpharep>=0),lines(rmax,gammarep,lwd=2,lty=2))
with(subset(thedat),lines(rmax,gammalo,lwd=2,lty=1))
dev.off()

# only if needs the tex output
# tikz('figs/figAnnex-rm-alpha.tex', width=3, height=3.5, standAlone=FALSE)
# with(thedat,plot(rmax,alphahi,type="l",lwd=2,ylim=c(0,3.5),xlim=c(1.2,2.4),xlab="",ylab="CoESS virulence ($\\alpha^*$)"))
# with(subset(thedat,alpharep>0),lines(rmax,alpharep,lwd=2,lty=2))
# with(subset(thedat),lines(rmax,alphalo,lwd=2,lty=1))
# dev.off()      
# 
# tikz('figs/figAnnex-rm-gamma.tex', width=3, height=3.5)
# with(thedat,plot(rmax,gammahi,type="l",lwd=2,ylim=c(0,4.5),xlab="Maximal reproduction rate ($r_{m}$)",ylab="CoESS recovery rate ($\\gamma^*$)"))
# with(subset(thedat,alpharep>=0),lines(rmax,gammarep,lwd=2,lty=2))
# with(subset(thedat),lines(rmax,gammalo,lwd=2,lty=1))
# dev.off()      