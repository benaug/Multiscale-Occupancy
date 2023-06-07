#this test script assumes you are using 3D data, i.e. p varies over replicates

#Data dimensions (definitions below may be specific to eDNA interpretation)
M <- 100 #sites
J <- 10 #samples
K <- 5 #replicates

#nimble code set up for unequal usage. Plugging in equal usage here
J1D <- rep(J,M)
K2D <- matrix(K,M,J)

#Here is an example of unequal usage
J1D <- sample(c(5,10,15),M,replace=TRUE)
K2D <- matrix(sample(c(3,5,8),M*max(J1D),replace=TRUE),M,max(J1D))


#parameter values
psi <- 0.25
theta <-  0.25
p <- 0.5

#simulate data
z <- rbinom(M,1,psi) #site occupancy states
w <- matrix(0,M,max(J1D)) #sample occupancy states
y <- array(0,dim=c(M,max(J1D),max(K2D))) #detections
for(i in 1:M){ #sites
  if(z[i]==1){
    for(j in 1:J1D[i]){ #samples
      w[i,j] <- rbinom(1,1,theta)
      y[i,j,1:K2D[i,j]] <- rbinom(K2D[i,j],1,p*w[i,j]) #replicates
    }
  }
}

#inits
w.init <- 1*(apply(y,c(1,2),sum)>0)
z.init <- 1*(rowSums(w.init)>0)

library(nimble)
library(coda)
source("State Sampler 3D.R")

#Nimble model
#note, to use "state sampler 2D", you *must* have objects "psi.site", "theta.site.sample", and "p.site.sample.rep"
#plugging in fixed parameter to these structures here, but can let them vary at these levels
NimModel <- nimbleCode({
  logit(psi) ~ dlogis(0,1) #P(measurable C at site)
  logit(theta) ~ dlogis(0,1) #P(measurable C in sample),
  logit(p) ~ dlogis(0,1) #P(measurable C in replicate)
  for(i in 1:M){
    psi.site[i] <- psi
    z[i] ~ dbern(psi) #is site occupied by edna?
    for(j in 1:J[i]){
      theta.site.sample[i,j] <- theta  #can multiply by usage here, e.g., theta*usage[i,j], where usage is 0 or 1
      w[i,j] ~ dbern(theta.site.sample[i,j]*z[i]) #is there edna in this sample?
      for(k in 1:K[i,j]){
        p.site.sample.rep[i,j,k] <- p #can multiply by usage here
        y[i,j,k] ~ dbinom(p=p.site.sample.rep[i,j,k]*w[i,j],size=1)
      }
    }
  }
})# end model


Niminits <- list(z=z.init,w=w.init)
constants <- list(M=M,J=J1D,K=K2D)
Nimdata <- list(y=y)

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
                type = 'StateSampler3D',control = list(M=M,J=J1D,K=K2D,z.det=z.det),
                silent = TRUE)
#nimble will tell you no samplers are assigned to w, because I set the target to z only. Ignore.

# Build and compile
Rmcmc <- buildMCMC(conf)
# runMCMC(Rmcmc,niter=10) #used for debugging/modifying custom sampler (put browser() in there)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


Cmcmc$run(5000,reset=FALSE)  #can keep running this line to continue run
mvSamples = as.matrix(Cmcmc$mvSamples)

burnin=500
plot(coda::mcmc(mvSamples[-c(1:burnin),]))


