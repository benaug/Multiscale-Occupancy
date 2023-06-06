## sampler to update z's and w's
StateSampler2D <- nimbleFunction(
  contains = sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # Defined stuff
    M <- control$M
    J <- control$J
    K <- control$K
    z.det <- control$z.det
    calcNodes <- model$getDependencies(c("z","w",'y')) #stuff to copy after update
  },
  run = function() {
    psi <- model$psi.site
    theta <- model$theta.site.sample
    p <- model$p.site.sample
    # z.det <- rep(0,M)
    z.preds <- rep(0,M)
    w.preds <- matrix(0,M,J)
    for(i in 1:M){
      # z.det[i] <- 1*(sum(model$y[i,])>0) #can supply via control, but calculating here so I can do parallel without recompiling
      lp.z <- rep(0,2)
      lp.w <- matrix(0,J,2)
      lp.w.total <- rep(0,J)
      if(z.det[i]==1){
        lp.z[1] <- -Inf #cannot have z=0 if there were detections
      }else{
        lp.z[1] <- dbinom(0,1,psi[i],log=TRUE)#transect level occupancy
      }
      lp.z[2] <- dbinom(1,1,psi[i],log=TRUE)
      #w stuff for only if z=1
      for(j in 1:J){
        for(w.state in 1:2){
          #w logProb
          lp.w[j,w.state] <- dbinom(w.state-1,1,theta[i,j],log=TRUE)
          #add y logProb
          lp.w[j,w.state] <- lp.w[j,w.state] + dbinom(model$y[i,j],K,p[i,j]*(w.state-1),log=TRUE)
        }
        max.lp <- max(lp.w[j,])
        lp.w.total[j] <- max.lp+log(sum(exp(lp.w[j,]-max.lp)))
      }
      #only need to add lp.w to lp.z[2]
      lp.z[2] <- lp.z[2]+sum(lp.w.total)
      if(z.det[i]>0){
        z.preds[i] <- 1
      }else{
        z.probs <- exp(lp.z)
        z.probs <- z.probs/sum(z.probs)
        z.preds[i] <- rcat(1,z.probs)-1
      }
      if(z.preds[i]==1){
        #sample w_ij's if z_i=1
        for(j in 1:J){
            w.probs <- exp(lp.w[j,])
            w.probs <- w.probs/sum(w.probs)
            w.preds[i,j] <- rcat(1,w.probs)-1
        }
      }
    }
    model$z <<- z.preds
    model$w <<- w.preds
    model$calculate(calcNodes)
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {} )
)
