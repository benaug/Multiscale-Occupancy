#this test script assumes you are using 2D data, i.e. p does not vary over replicates
#so y is summed over replicates

#Data dimensions (definitions below may be specific to eDNA interpretation)
M <- 100 #sites
J <- 10 #samples
K <- 5 #replicates

#parameter values
psi <- 0.25
theta <-  0.25
p <- 0.5

#simulate data
z <- rbinom(M,1,psi) #site occupancy states
w <- matrix(0,M,J) #sample occupancy states
y <- array(0,dim=c(M,J,K)) #detections
for(i in 1:M){ #sites
  if(z[i]==1){
    for(j in 1:J){ #samples
      w[i,j] <- rbinom(1,1,theta)
      y[i,j,] <- rbinom(K,1,p*w[i,j]) #replicates
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
#note, to use "state sampler 2D", you *must* have objects "psi.site", "theta.site", and "p.site.sample"
#plugging in fixed parameter to these structures here, but can let them vary at these levels
NimModel <- nimbleCode({
  logit(psi) ~ dlogis(0,1) #P(measurable C at site)
  logit(theta) ~ dlogis(0,1) #P(measurable C in sample),
  logit(p) ~ dlogis(0,1) #P(measurable C in replicate)
  for(i in 1:M){
    psi.site[i] <- psi
    theta.site[i] <- theta
    z[i] ~ dbern(psi) #is site occupied by edna?
    for(j in 1:J){
      p.site.sample[i,j] <- p
      w[i,j] ~ dbern(theta.site[i]*z[i]) #is there edna in this sample?
      y[i,j] ~ dbinom(p=p.site.sample[i,j]*w[i,j],size=K)
    }
  }
})# end model


Niminits <- list(z=z.init,w=w.init)
constants <- list(M=M,J=J,K=K)
Nimdata <- list(y=y2D)

# set parameters to monitor
parameters <- c('psi','theta','p')

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
                type = 'StateSampler2D',control = list(M=M,J=J,K=K,z.det=z.det),
                silent = TRUE)
#nimble will tell you no samplers are assigned to w, because I set the target to z only. Ignore.

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #used for debugging/modifying custom sampler (put browser() in there)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


Cmcmc$run(5000,reset=TRUE)
mvSamples = as.matrix(Cmcmc$mvSamples)

burnin=500
plot(coda::mcmc(mvSamples[-c(1:burnin),]))


