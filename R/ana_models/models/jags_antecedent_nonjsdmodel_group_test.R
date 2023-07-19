model{
    for(k in 1:n.species){  # Loop through species
      for(i in 1:n.sites){ #loop through the sites
      # Likelihood
        # Data are around occupancy probability
        #this is "apparent" occupancy because it's not
        # corrected for detection
        y[k,i] ~ dbern(p[k,i])
        
        #occupancy dependent on covariates
        #logistic model of occupancy
        logit(p[k,i]) <- a0[GroupID[k]] + 
          a1[GroupID[k]]*SoilMAnt[GroupID[k],i] 
        
        
        #could add in random effects and other covariates to above
        # if desired
        
      } #n.sites

    } #n.species
  
  #Species-group level Priors
  for(g in 1:n.groups){
    a0[g] ~ dnorm(mu.a0, tau.a0)
    a1[g] ~ dnorm(mu.smant, tau.smant)
  }
  
  #Antecedent climate values
  for(g in 1:n.groups){
  for(i in 1:n.sites){
    #antecedent climate value is the sum of all the
    #weighted lagged values
    SoilMAnt[g, i] <- sum(SoilMTemp[g,i,])

    #these are where the data you input into the model
    #come in for soil moisutre - here they're multiplied
    #by the weights the model has given them for their
    #relative importance in driving the climate effect
    for(t in 1:n.lag){
      SoilMTemp[g,i,t] <- SoilM[i,t]*wA[g,t]
    }
  }
}

  #Employing "delta trick" to give vector of weights dirichlet priors
  #this is doing the dirichlet in two steps
  #see Ogle et al. 2015 SAM model paper in Ecology Letters
  for(g in 1:n.groups){
    #Antecedent climate priors
    #sum of the weights for the climate lag, so we
    # can divide the total value by 1
    sumA[g] <- sum(deltaA[g,])
    
  for(t in 1:n.lag){ #for total number of lags
   #the weigths for soilmoisture - getting them to sum to 1
    wA[g,t] <- deltaA[g,t]/sumA[g]
    #both follow relatively uninformative gamma priors
    deltaA[g,t] ~ dgamma(1,1)
  }
  }
  
  # Hyperpriors (community level)
  #easier to make uniform (~beta(1,1)) prior on non-logistic scale
  a0.mean ~ dbeta(1,1)
  mu.a0 <- logit(a0.mean)
  sd.a0 ~ dunif(0, 5)
  tau.a0 <- pow(sd.a0, -2)

  mu.smant ~ dunif(-5,5)
  sd.smant ~ dunif(0,5)
  tau.smant <- pow(sd.smant,-2)

  
  }
  