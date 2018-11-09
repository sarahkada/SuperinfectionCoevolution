rm(list=ls(all=TRUE)) 
source("supercoev-superfunction.R")
require(tikzDevice)
#############################################################################
############## alpha,gamma vs sigma' (called sigmap here) - fig 7b and 7d

# drawing only:
thedat <- read.table(file = "dat/fig7-vsigmaprime.csv", header=T)

epsoutput("figs/fig7b-vsigmaprim-alpha.eps")
with(thedat,plot(sigmap,alphahi,type="l",lwd=2,xlim=c(-0.02,0.4), ylim=c(0,7.5),xlab=expression("Strengh of within-host competition" ~ (sigma^prime)),ylab=expression("CoESS virulence" ~ (alpha)), cex.lab=1.1))
with(subset(thedat,alpharep>0),lines(sigmap,alpharep,lwd=2,lty=2))
with(subset(thedat),lines(sigmap,alphalo,lwd=2,lty=1))
dev.off()

epsoutput("figs/fig7b-vsigmaprim-gamma.eps")
with(thedat,plot(sigmap,gammahi,type="l",lwd=2,ylim=c(0,4.5),xlab=expression("Strengh of within-host competition" ~ (sigma^prime)),ylab=expression("CoESS recovery rate" ~ (gamma)), cex.lab=1.1))
with(subset(thedat,gammarep>0),lines(sigmap,gammarep,lwd=2,lty=2))
with(subset(thedat),lines(sigmap,gammalo,lwd=2,lty=1))
dev.off()


#### or re-run the simulations:
kappa<-0.01
betam<-10
mu<-1
rmax<-2
epsilon<-1
c<-0.05
smax<-10
cup<-0.0
sigma<- 0.01

# sigmap from 0 to -0.2 and high coess and reppeller meet
thedat2<-NULL
sigmap<- 0.011
parms<-c(sigma=sigma, sigmap=sigmap, smax=smax,kappa=kappa, betam=betam, mu=mu, rmax=rmax,epsilon=epsilon, c=c, cup=cup)
hicoess  <- multiroot(f = selgrad, start = c(3,4),positive=T)
eqhi<-equilibrium(c((hicoess$root)[1],(hicoess$root)[2]))
repellor <- multiroot(f = selgrad, start = c(1,0),positive=T)
locoess  <- uniroot(selgraduni, c(1,2), tol = 0.0001, maxiter = 1000)
eqlo<-equilibrium(c(locoess$root,0))

#thedat2<-rbind(thedat2,data.frame(sigmap=sigmap,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=(repellor$root)[1],gammarep=(repellor$root[2]),alphalo=locoess$root,gammalo=0,eqhiS=eqhi[1],eqhiI=eqhi[2],eqloS=eqlo[1],eqloI=eqlo[2]))

computeall=T

for (sigmap in c(seq(0.011, -0.025, by= - 0.0005), seq(-0.025, -0.0273, by= - 0.0001))) {  
  
  parms<-c(sigma=sigma, sigmap=sigmap, smax=smax,kappa=kappa, betam=betam, mu=mu, rmax=rmax,epsilon=epsilon, c=c, cup=cup)
  
  # locoess  <- uniroot(selgraduni, c(locoess$root-0.5,locoess$root+0.5), tol = 0.0001, maxiter = 1000)
  
  if(computeall) {
    repellor <- multiroot(f = selgrad, start = c((repellor$root)[1],(repellor$root)[2]),positive=T)
    hicoess  <- multiroot(f = selgrad, start = c((hicoess$root)[1],(hicoess$root)[2]),positive=T)
  }
  #  } else {
  #    repellor$root=c(-5,NA)
  #     hicoess$root=c(-1,NA)
  #   }
  if(abs((repellor$root)[1]-(hicoess$root)[1])<0.05)
  {
    computeall=F
  }
  
  thisdat <- data.frame(sigmap=sigmap,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=(repellor$root)[1],gammarep=(repellor$root[2]),alphalo=locoess$root,gammalo=0,eqhiS=eqhi[1],eqhiI=eqhi[2],eqloS=eqlo[1],eqloI=eqlo[2])
  
  thedat2<-rbind(thedat2,thisdat)
  
  #write.table(thisdat,paste("fig7-vssigmaprime",".dat",sep=""),append=T,row.names=F,col.names=F)
  
}

thedat<-NULL
sigmap<-  0.011
parms<-c(sigma=sigma, sigmap=sigmap, smax=smax,kappa=kappa, betam=betam, mu=mu, rmax=rmax,epsilon=epsilon, c=c, cup=cup)
hicoess  <- multiroot(f = selgrad, start = c(2,3),positive=T)
eqhi<-equilibrium(c((hicoess$root)[1],(hicoess$root)[2]))
repellor <- multiroot(f = selgrad, start = c(1,0),positive=T)
locoess  <- uniroot(selgraduni, c(1,2), tol = 0.0001, maxiter = 1000)
eqlo<-equilibrium(c(locoess$root,0))

#thedat<-rbind(thedat,data.frame(sigmap=sigmap,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=(repellor$root)[1],gammarep=(repellor$root[2]),alphalo=locoess$root,gammalo=0,eqhiS=eqhi[1],eqhiI=eqhi[2],eqloS=eqlo[1],eqloI=eqlo[2]))

computeall=T

for (sigmap in c(seq(0.012, 0.15,by=0.001), seq(0.15, 0.3,by=0.01), seq(0.3, 0.35, by=0.001))) {  
  
  parms<-c(sigma=sigma, sigmap=sigmap, smax=smax,kappa=kappa, betam=betam, mu=mu, rmax=rmax,epsilon=epsilon, c=c, cup=cup)
  
  hicoess  <- multiroot(f = selgrad, start = c((hicoess$root)[1],(hicoess$root)[2]),positive=T)
  eqhi<-equilibrium(c((hicoess$root)[1],(hicoess$root)[2]))
  
  # if(computeall) {
  #   repellor <- multiroot(f = selgrad, start = c((repellor$root)[1],(repellor$root)[2]),positive=T)
  #   locoess  <- uniroot(selgraduni, c(locoess$root-0.5,locoess$root+0.5), tol = 0.0001, maxiter = 1000)
  #   eqlo<-equilibrium(c(locoess$root,0))
  # }
  # else {
  #   repellor$root=c(-5,NA)
  #   locoess$root=-1
  #   eqlo<-c(NA,NA)
  # }
  
  
  # if(abs((repellor$root)[1]-(locoess$root))<0.01)
  # {
  #   computeall=F
  # }
  
  thisdat <- data.frame(sigmap=sigmap,alphahi=(hicoess$root)[1],gammahi=(hicoess$root)[2],alpharep=NA,gammarep=NA,alphalo=NA,gammalo=NA,eqhiS=eqhi[1],eqhiI=eqhi[2],eqloS=NA,eqloI=NA)
  
  thedat<-rbind(thedat,thisdat)
  
  #write.table(thisdat,paste("fig7-vssigmaprime",".dat",sep=""),append=T,row.names=F,col.names=F)
  
}

thedat2<-thedat2[order(thedat2$sigmap),]
thedat<-rbind(thedat2,thedat) 

#write.table(thedat,"dat/fig7-vsigmaprime.csv")

epsoutput("figs/fig7b-vsigmaprim-alpha.eps")
with(thedat,plot(sigmap,alphahi,type="l",lwd=2,xlim=c(-0.02,0.4), ylim=c(0,7.5),xlab=expression("Strengh of within-host competition" ~ (sigma^prime)),ylab=expression("CoESS virulence" ~ (alpha)), cex.lab=1.1))
with(subset(thedat,alpharep>0),lines(sigmap,alpharep,lwd=2,lty=2))
with(subset(thedat),lines(sigmap,alphalo,lwd=2,lty=1))
dev.off()

epsoutput("figs/fig7b-vsigmaprim-gamma.eps")
with(thedat,plot(sigmap,gammahi,type="l",lwd=2,ylim=c(0,4.5),xlab=expression("Strengh of within-host competition" ~ (sigma^prime)),ylab=expression("CoESS recovery rate" ~ (gamma)), cex.lab=1.1))
with(subset(thedat,gammarep>0),lines(sigmap,gammarep,lwd=2,lty=2))
with(subset(thedat),lines(sigmap,gammalo,lwd=2,lty=1))
dev.off()

# draw figure in tex format:
# tikz('figs/fig7b-vsigmaprim-alpha.tex', width=3, height=3.5)
# with(thedat,plot(sigmap,alphahi,type="l",lwd=2, c(-0.05,0.35), ylim=c(0,7.5),xlab="Within-host competition strengh ($\\sigma'$)",ylab="CoESS virulence ($\\alpha^*$)"))
# with(subset(thedat,alpharep>0),lines(sigmap,alpharep,lwd=2,lty=2))
# with(subset(thedat),lines(sigmap,alphalo,lwd=2,lty=1))
# dev.off()            
# 
# tikz('figs/fig7b-vsigmaprim-gamma.tex', width=3, height=3.5)
# with(thedat,plot(sigmap,gammahi,type="l",lwd=2,xlim=c(-0.05,0.35),ylim=c(0,4.5),xlab="Within-host competition strengh ($\\sigma'$)",ylab="CoESS recovery rate ($\\gamma^*$)"))
# with(subset(thedat,gammarep>0),lines(sigmap,gammarep,lwd=2,lty=2))
# with(subset(thedat),lines(sigmap,gammalo,lwd=2,lty=1))
# dev.off()      