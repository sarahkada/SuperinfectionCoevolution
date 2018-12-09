rm(list=ls(all=TRUE)) # removes all stored characters
cat("\014") # clear console

source("supercoev-isoclines.R")

kappa<-0.01
betam<-10
mu<-1
rmax<-2
sigma0<-0.0
epsilon<-1
c<-0.05
cup<-0.0
smtg<- 1
parms<-c(smtg=smtg, kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
thedat<-NULL

# the idea here is to find the host and parasite isocline: finding the solution for alpha when gamma varies and conversly. 
# this numerical computation is very sensitive, so we have to proceed step-wise 
#########################################################
#    1. Parasite isocline, with varying, gamma sigma=0  #
#########################################################
# Failed with g_step = 3.1
psus <- 0
pinf <- 0
a = c()
b = c()

for (g_step in seq(0.0,4.0,by=0.1)){
  #  "tryCatch" makes the iteration continue even is uniroot failed
  # Example of uniroot trouble
  #		https://stat.ethz.ch/pipermail/r-help/2007-February/124692.html
  #		http://www.student.tue.nl/V/j.g.v.d.pol/Teaching/R/tutorial.asp (tout a la fin)
  #	Many things can go wrong when implementing the function uniroot. Perhaps the error that 
  #	occurs most often is "Error: Values at end points not of opposite sign". This has to do 
  #	with how uniroot is implemented: it expects to be fed a function that is negative on 
  #	one side of the zero, and positive on the other.
  tryCatch(
{
  iso <- uniroot(selgradparalpha, interval=c(0.0,2.0), gam=g_step)
  
  a <- c(a, iso$root)
  b <- c(b, g_step)
  
  thedat <- rbind(thedat, data.frame(alphaiso=a, gamma=b))
},
error=function(cond) {
  message(cond)
  message(cat(" | Failed with g_step = ", g_step, "\n", sep=""))
},
warning=function(cond) {
  message(cond)
},
finally={
}
  )    
}
#Fails with g_step = 5.3 so re-run for gamma = 4 to 8
aa = c()
bb = c()
thedat2<-NULL

for (g_step in seq(3.0,6.0,by=0.1)){
  tryCatch(
{
  iso <- uniroot(selgradparalpha, interval=c(1.5,2.5), gam=g_step)
  
  aa <- c(aa, iso$root)
  bb <- c(bb, g_step)
  
  thedat2 <- rbind(thedat2, data.frame(alphaiso=aa, gamma=bb))
},
error=function(cond) {
  message(cond)
  message(cat(" | Failed with g_step = ", g_step, "\n", sep=""))
},
warning=function(cond) {
  message(cond)
},
finally={
}
  )    
}
a <- c(a,aa)
b <- c(b,bb)

#########################################################
#     1.2           step on alpha with sigma=1         #
#########################################################
####### Parameters
sigma0<-1.0
parms<-c(smtg=smtg, kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
thedat2<-NULL

psus <- 0
pinf <- 0
d = c()
e = c()

for (g_step in seq(0.0,4.5,by=0.1)) {
  tryCatch(
{
  iso <- uniroot(selgradparalpha, interval=c(1.5,3.0), gam=g_step)#, tol = 0.0001, maxiter = 1000)
  
  d <- c(d, iso$root)
  e <- c(e, g_step)
  thedat2 <- rbind(thedat2, data.frame(alphaiso=d, gamma=e))
},
error=function(cond) {
  message(cond)
  message(cat(" | Failed with g_step = ", g_step, "\n", sep=""))
},
warning=function(cond) {
  message(cond)
},
finally={
}
  )    
}

# loop for gamma from 4 to 8 (fails at g_step= 6)
dd = c()
ee = c()

for (g_step in seq(4.5,8.0,by=0.1)) {
  tryCatch(
{
  iso <- uniroot(selgradparalpha, interval=c(2.0,3.0), gam=g_step)#, tol = 0.0001, maxiter = 1000)
  
  dd <- c(dd, iso$root)
  ee <- c(ee, g_step)
},
error=function(cond) {
  message(cond)
  message(cat(" | Failed with g_step = ", g_step, "\n", sep=""))
},
warning=function(cond) {
  message(cond)
},
finally={
}
  )    
}
d <- c(d,dd)
e <- c(e,ee)

#########################################################
#   1.3         Step on alpha, sigma=6                ##
#########################################################
####### Parameters
sigma0<-6.0
parms<-c(smtg=smtg, kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
thedat3<-NULL

psus <- 0
pinf <- 0
f = c()
g = c()

for (g_step in seq(0.0,8.0,by=0.1)) {
  tryCatch(
{
  iso <- uniroot(selgradparalpha, interval=c(3.0,4.0), gam=g_step)#, tol = 0.0001, maxiter = 1000)
  
  f <- c(f, iso$root)
  g <- c(g, g_step)
  
  thedat3 <- rbind(thedat3, data.frame(alphaiso=f, gamma=g))
},
error=function(cond) {
  message(cond)
  message(cat(" | Failed with g_step = ", g_step, "\n", sep=""))
},
warning=function(cond) {
  message(cond)
},
finally={
}
  )    
} 

############### run for gamma > 4
thedat3<-NULL
psus <- 0
pinf <- 0
ff = c()
gg = c()

for (g_step in seq(4.8,8.0,by=0.1)) {
  tryCatch(
{
  iso <- uniroot(selgradparalpha, interval=c(3.0,5.0), gam=g_step)#, tol = 0.0001, maxiter = 1000)
  
  ff <- c(ff, iso$root)
  gg <- c(gg, g_step)
  
  thedat3 <- rbind(thedat3, data.frame(alphaiso=ff, gamma=gg))
},
error=function(cond) {
  message(cond)
  message(cat(" | Failed with g_step = ", g_step, "\n", sep=""))
},
warning=function(cond) {
  message(cond)
},
finally={
}
  )    
}

f <- c(f,ff)
g <- c(g,gg)

#########################################################
#   2.   Host isocline with varying alpha             #
#########################################################
sigma0<-0.0
parms<-c(smtg=smtg, kappa=kappa, betam=betam, mu=mu, rmax=rmax, sigma0=sigma0, epsilon=epsilon, c=c, cup=cup)
thedat<-NULL

psus <- 0
pinf <- 0
ab = c()
ba = c()
thedat = data.frame()
#  fails at a_step= 2.95 
for (a_step in seq(0.05,8.0,by=0.01)) {
  tryCatch(
{
  iso <- uniroot(selgradhostalpha, interval=c(0.0,12.0), alp=a_step)
  
  ab <- c(ab, iso$root)
  ba <- c(ba, a_step)
  
  thedat <- rbind(thedat, data.frame(alpha=ba, gammaiso=ab))
},
error=function(cond) {
  message(cond)
  message(cat(" | Failed with a_step = ", a_step, "\n", sep=""))
},
warning=function(cond) {
  message(cond)
},
finally={
}
  )    
}

# partie de alpha=7 à 8
psus <- 0
pinf <- 0
ab_end = c()
ba_end = c()
iso <- uniroot(selgradhostalpha, interval=c(0,1), alp=7.04)#, tol = 0.0001, maxiter = 1000)
for (a_step in seq(7.5,7.7,by=0.05)) {
  tryCatch(
{
  iso <- uniroot(selgradhostalpha, interval=c(iso$root-0.5,iso$root+0.5), alp=a_step)#, tol = 0.0001, maxiter = 1000)
  
  ab_end <- c(ab_end, iso$root)
  ba_end <- c(ba_end, a_step)
  
},
error=function(cond) {
  message(cond)
  message(cat(" | Failed with a_step = ", a_step, "\n", sep=""))
},
warning=function(cond) {
  message(cond)
},
finally={
}
  )    
}

ab <- c(ab,ab_end)
ba <- c(ba,ba_end)

#### alpha= 1 to 1.31
s = c()
t = c()

for (a_step in seq(1.0,1.31,by=0.01)) { 
  tryCatch(
{
  iso <- uniroot(selgradhostalpha, interval=c(0.0,1.0), alp=a_step)
  
  s <- c(s, iso$root)
  t <- c(t, a_step)
  
  thedat <- rbind(thedat, data.frame(alpha=s, gammaiso=t))
},
error=function(cond) {
  message(cond)
  message(cat(" | Failed with a_step = ", a_step, "\n", sep=""))
},
warning=function(cond) {
  message(cond)
},
finally={
}
  )    
}

#fails at the beginning 
y = c()
z = c()
for (a_step in seq(1.31,1.36,by=0.01)) {
  tryCatch(
{
  iso <- uniroot(selgradhostalpha, interval=c(0.0,2.0), alp=a_step)#, tol = 0.0001, maxiter = 1000)

  y <- c(y, iso$root)
  z <- c(z, a_step)
},
error=function(cond) {
  message(cond)
  message(cat(" | Failed with a_step = ", a_step, "\n", sep=""))
},
warning=function(cond) {
  message(cond)
},
finally={
}
  )    
}

z <- c(t,z)
y <- c(s,y)
#########################################################
#                   plot figure                      ###
#########################################################
library(shape)
epsoutput("figs/fig1-isoclines.eps")
# require(tikzDevice)
# tikz('figs/fig1-isoclines.tex', width=4, height=3.5, standAlone=T)# 
plot(a,b,type="l",lwd=2,xlim=c(0,8),ylim=c(0,5),xlab="ES virulence ($\\alpha$)",ylab="ES recovery rate ($\\gamma$)") # parasite isocline sigma=0
lines(d,e,lwd=2,lty=2) # parasite isocline, sigma=1
lines(f,g,type="l",lwd=2, lty=3) # parasite isocline, sigma=6
lines(ba,ab,type="l",lwd=3) # host isocline
lines(z,y,type="l",lwd=3) # host isocline for low alpha 
points(2.1,3.4, pch = 20)
points(2.65,3.52, pch = 20)
points(3.85,3.08, pch = 20)
points(1,0, pch = 20)
points(1.27,0.6, pch = 21, bg = "white", lwd=3,  cex = 0.8)
text(3.3, 4.0, paste("$\\sigma = 1$"), cex = 0.7)
text(4.8, 4.0, paste("$\\sigma = 6$"), cex = 0.7)
text(1.5, 4.0, paste("$\\sigma  = 0$"), cex = 0.7)
lines(c(0, 1.25), c(0,0),lwd=3)
lines(c(7.722, 8.5), c(0,0),lwd=3)
dev.off()
