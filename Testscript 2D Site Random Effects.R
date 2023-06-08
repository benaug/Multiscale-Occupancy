#this test script assumes you are using 2D data, i.e. p does not vary over replicates
#so y is summed over replicates. See 3D test script if p varies across replicates.
#this script is also set up for site-level random effects for theta and p
#uses non-centered parameterization

#Data dimensions (definitions below may be specific to eDNA interpretation)
M <- 100 #sites
J <- 10 #samples
K <- 10 #replicates

#nimble code set up for unequal usage. Plugging in equal usage here
J1D <- rep(J,M)
K2D <- matrix(K,M,J)

# #Here is an example of unequal usage
# J1D <- sample(c(5,10,15),M,replace=TRUE)
# K2D <- matrix(sample(c(3,5,8),M*max(J1D),replace=TRUE),M,max(J1D))


#parameter values
psi <- 0.25
theta <-  0.25
p <- 0.25
sd.theta <- 0.5 #theta sd on logit scale
sd.p <- 0.5 #p sd on logit scale
#simulate data
z <- rbinom(M,1,psi) #site occupancy states
w <- matrix(0,M,max(J1D)) #sample occupancy states
y <- array(0,dim=c(M,max(J1D),max(K2D))) #detections
for(i in 1:M){ #sites
  if(z[i]==1){
    theta.site=plogis(rnorm(1,qlogis(theta),sd.theta))
    p.site=plogis(rnorm(1,qlogis(p),sd.p))
    for(j in 1:J1D[i]){ #samples
      w[i,j] <- rbinom(1,1,theta.site)
      y[i,j,1:K2D[i,j]] <- rbinom(K2D[i,j],1,p.site*w[i,j]) #replicates
    }
  }
}

#data summed over replicates
y2D <- apply(y,c(1,2),sum)

#inits
w.init <- 1*(apply(y,c(1,2),sum)>0)
z.init <- 1*(rowSums(w.init)>0)

library(nimble)
library(coda)
source("State Sampler 2D.R")

#Nimble model
#note, to use "state sampler 2D", you *must* have objects "psi.site", "theta.site.sample", and "p.site.sample"
#plugging in fixed parameter to these structures here, but can let them vary at these levels
NimModel <- nimbleCode({
  logit(psi) ~ dlogis(0,1) #P(measurable C at site)
  logit(theta) ~ dlogis(0,1) #P(measurable C in sample),
  logit(p) ~ dlogis(0,1) #P(measurable C in replicate)
  theta.sd ~ T(dt(0, df = 7, sigma = 1),0,) #half-t prior
  p.sd ~ T(dt(0, df = 7, sigma = 1),0,) #half-t prior
  # theta.sd ~ dunif(0,5) #uniform prior
  # p.sd ~ dunif(0,5) #unifor prior
  for(i in 1:M){
    psi.site[i] <- psi
    z[i] ~ dbern(psi) #is site occupied by edna?
    #non-centered random effects for theta and p
    theta.unit[i] ~ dnorm(0,sd=1)
    logit(theta.site[i]) <- logit(theta) + theta.unit[i]*theta.sd
    p.unit[i] ~ dnorm(0,sd=1)
    logit(p.site[i]) <- logit(p) + p.unit[i]*p.sd
    for(j in 1:J[i]){
      theta.site.sample[i,j] <- theta.site[i]
      p.site.sample[i,j] <- p.site[i]
      w[i,j] ~ dbern(theta.site.sample[i,j]*z[i]) #is there edna in this sample?
      y[i,j] ~ dbinom(p=p.site.sample[i,j]*w[i,j],size=K[i,j])
    }
  }
})# end model

Niminits <- list(z=z.init,w=w.init)
constants <- list(M=M,J=J1D,K=K2D)
Nimdata <- list(y=y2D)

# set parameters to monitor
parameters <- c('psi','theta','p','theta.sd','p.sd')

Rmodel <- nimbleModel(code=NimModel, constants=constants, data=Nimdata,check=FALSE,
                      inits=Niminits)
conf <- configureMCMC(Rmodel,monitors=parameters, thin=1,useConjugacy = FALSE)

#add custom w and z sampler
#can skip this replacement to run model with normal binary samplers assigned by nimble
#precalculate z.det - which sites have a detection
w.det <- 1*(apply(y,c(1,2),sum)>0)
z.det <- 1*(rowSums(w.det)>0)
conf$removeSampler(c('z','w'))
conf$addSampler(target = paste("z[1:",M,"]"),
                type = 'StateSampler2D',control = list(M=M,J=J1D,K=K2D,z.det=z.det),
                silent = TRUE)
#nimble will tell you no samplers are assigned to w, because I set the target to z only. Ignore.

#might want to block these--posteriors can be strongly correlated
# conf$removeSampler(c('logit_theta','theta.sd'))
# conf$addSampler(target = c('logit_theta','theta.sd'),type = 'RW_block',control = list(),silent = TRUE)
# conf$removeSampler(c('logit_p','p.sd'))
# conf$addSampler(target = c('logit_p','p.sd'),type = 'RW_block',control = list(),silent = TRUE)


# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #used for debugging/modifying custom sampler (put browser() in there)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


Cmcmc$run(15000,reset=FALSE) #can keep running this line to continue run
mvSamples = as.matrix(Cmcmc$mvSamples)

burnin=500
plot(coda::mcmc(mvSamples[-c(1:burnin),]))
