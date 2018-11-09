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

# epidemiological model
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

# selection gradient for the parasite, depends on alpha and gamma
selgradpar <- function(alpha,gamma) {
	   dbeta(alpha)-beta(alpha)/(mu+alpha+gamma+sigma(alpha,gamma)*beta(alpha)*pinf)
}

# selection gradient for the host, depends on alpha and gamma
selgradhost <- function(alpha,gamma) {
	   #expression sans densité dépendance # drS(gamma)+drI(gamma)*beta(alpha)*pinf/(mu+alpha+gamma)+(rS(gamma)-mu)/(mu+alpha+gamma)
  drS(gamma)*(1-kappa*(psus+pinf))+drI(gamma)*(1-kappa*(psus+pinf))*beta(alpha)*pinf/(mu+alpha+gamma)+(rS(gamma)*(1-kappa*(psus+pinf))-mu)/(mu+alpha+gamma)
}

selgrad <- function(x) {
  cat("alpha=", x[1], " | gamma=", x[2], "\n", sep="")
	# initialise solver
	yini <- c(pS=0.5,pI=0.001)
	pars <- c(parms, alpha=x[1], gamma = x[2])
	times <- seq(0,100,1)

	# integrate epidemiological model

	sol <- ode(y=yini, times=times, func=model, parms=pars, method="ode45")		
	
	psus <<- sol[,"pS"][101]
	pinf <<- sol[,"pI"][101]

	return(c(F1 = selgradhost(x[1],x[2]), F2 = selgradpar(x[1],x[2])))

#with(as.list(pars),c(alpha,gamma))
}

equilibrium <-function(x) {
	# initialise solver
	yini <- c(pS=0.5,pI=0.001)
	pars <- c(parms, alpha=x[1], gamma = x[2])
	times <- seq(0,100,1)

	# integrate epidemiological model

	sol <- ode(y=yini, times=times, func=model, parms=pars,method="ode45")		
	
	psus <<- sol[,"pS"][101]
	pinf <<- sol[,"pI"][101]

	c(psus,pinf)
} # returns the equilibrium number of infected and susceptibles

selgraduni <- function(x) {
	# initialise solver
	yini <- c(pS=0.5,pI=0.001)
	pars <- c(parms, alpha=x[1], gamma = 0)
	times <- seq(0,100,1)
	# integrate epidemiological model

	sol <- ode(y=yini, times=times, func=model, parms=pars, method="ode45")		
	
	psus <<- sol[,"pS"][101]
	pinf <<- sol[,"pI"][101]
#return
	selgradpar(x,0)
}

epsoutput <- function(file) {
  setEPS()
  postscript(file)
}
