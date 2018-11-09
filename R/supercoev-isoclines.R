rm(list=ls(all=TRUE)) # removes all stored characters
library(deSolve)
library(rootSolve)

rS <- function(gamma) {
  rmax/(1+c*gamma)
}

rI <- function(gamma) {
  epsilon*rS(gamma)-cup*gamma
}

beta <- function(alpha) {
  betam*alpha/(1+alpha)
}

dbeta <- function(alpha) {
  betam/((1+alpha)*(1+alpha))
}

drS <- function(gamma) {
  -rmax*c/((1+c*gamma)*(1+c*gamma))
}

drI <- function(gamma) {
  epsilon*drS(gamma)-cup
}


model <- function(time, y, parms){
  with(as.list(c(y, parms)), {
    dpS <- (rS(gamma)*pS+rI(gamma)*pI)*(1-kappa*(pS+pI)) -mu*pS - beta(alpha)*pS*pI+gamma*pI
    dpI <- beta(alpha)*pS*pI-(mu+alpha+gamma)*pI
    list(c(dpS, dpI))
  })
}

sigma <- function(alpha,gamma) {
  sigma0
}

selgradpar <- function(alpha,gamma) {
  dbeta(alpha)-beta(alpha)/(mu+alpha+gamma+sigma(alpha,gamma)*beta(alpha)*pinf)
}

selgradhost <- function(alpha,gamma) {
 # drS(gamma)+drI(gamma)*beta(alpha)*pinf/(mu+alpha+gamma)+(rS(gamma)-mu)/(mu+alpha+gamma)
  drS(gamma)*(1-kappa*(psus+pinf))+drI(gamma)*(1-kappa*(psus+pinf))*beta(alpha)*pinf/(mu+alpha+gamma)+(rS(gamma)*(1-kappa*(psus+pinf))-mu)/(mu+alpha+gamma)
}

selgraduni <- function(x) {
  # initialise solver
  yini <- c(pS=0.5,pI=0.001)
  pars <- c(parms, alpha=x[1], gamma = 0)
  times <- seq(0,100,1)
 
  # integrate epidemiological model
  
  sol <- ode(y=yini, times=times, func=model, parms=pars)		
  
  psus <<- sol[,"pS"][101]
  pinf <<- sol[,"pI"][101]
  #return
  selgradpar(x,0)
}

selgradparalpha <- function(x, gam) {

  cat("gamma=", gam, " | x=", x, "\n", sep="")

  # initialise solver
  yini <- c(pS=0.5,pI=0.001)
  pars <- c(parms, alpha=x, gamma = gam)
  times <- seq(0,100,1)
  
  # integrate epidemiological model
  
  sol <- ode(y=yini, times=times, func=model, parms=pars, method="ode45")# préciser la méthode permet d'éviter les erreurs+faster
  
  psus <<- sol[,"pS"][101]
  pinf <<- sol[,"pI"][101]
  cat("pinf=", pinf, "\n", sep="")

  #return
  selgradpar(x, gam)
}

selgradhostalpha <- function(x, alp) {
  
  cat("alpha=", alp, " | x=", x, "\n", sep="")
  
  # initialise solver
  yini <- c(pS=0.5,pI=0.001)
  pars <- c(parms, alpha=alp, gamma = x)
  times <- seq(0,100,1)
  
  # integrate epidemiological model
  
  sol <- ode(y=yini, times=times, func=model, parms=pars, method="ode45")
  
  psus <- sol[,"pS"][101]
  pinf <<- sol[,"pI"][101]
  cat("pinf=", pinf, "\n", sep="")
  
  #return
  selgradhost(alp, x)
}

selgradhostgamma <- function(x, gam) {
  
  cat("gamma=", gam, " | x=", x, "\n", sep="")
  
  # initialise solver
  yini <- c(pS=0.5,pI=0.001)
  pars <- c(parms, alpha=x, gamma = gam)
  times <- seq(0,100,1)
  
  # integrate epidemiological model
  
  sol <- ode(y=yini, times=times, func=model, parms=pars, method="ode45")
  
  psus <- sol[,"pS"][101]
  pinf <<- sol[,"pI"][101]
  cat("pinf=", pinf, "\n", sep="")
  cat("psus=", psus, "\n", sep="")
  
  #return
  selgradhost(x, gam)
}

epsoutput <- function(file) {
  setEPS()
  postscript(file)
}