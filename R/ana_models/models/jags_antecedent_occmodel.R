model{
    for(k in 1:n.species){  # Loop through species
      for(i in 1:n.sites){ #loop through the sites
      # Likelihood
        # Data are around occupancy probability
        #this is "apparent" occupancy because it's not
        # corrected for detection
        y[k,i] ~ dbern(psi[k,i])
        
        #occupancy dependent on covariates
        logit(psi[k,i]) <- a0[k] + 
                            a1*TAnt[i] + 
                            a2*PAnt[i]
        
        #could add in random effects and other covariates to above
        # if desired
        #could also make the a1, a2 or other parameters dependent on
        # a speciesID or species group ID (e.g. effect of temperature
        # on grasses vs. forbes, large vs. small seeds, etc)
        
      } #n.sites

      
  #Species-level Priors
  a0[k] ~ dnorm(mu.a0, tau.a0)
  a1[k] ~ dnorm(mu.tant, tau.tant)
  a2[k] ~ dnorm(mu.ppt, tau.ppt)
  
    } #n.species
  
  #Antecedent climate values
  for(i in 1:n.sites){
    #antecedent climate value is the sum of all the 
    #weighted lagged values
    TAnt[i] <- sum(TmaxTemp[i,])
    Pant[i] <- sum(PPTTemp[i,])
    
    #these are where the data you input into the model
    #come in for temp and ppt - here they're multiplied
    #by the weights the model has given them for their 
    #relative importance in driving the climate effect
    for(t in 1:n.lag){
      TmaxTemp[i,t] <- Tmax[i,t]*wA[t]
      PPTTemp[i,t] <- PPT[i,t]*wB[t]
    }
    
  }

  #Antecedent climate priors
  #sum of the weights for the climate lag, so we 
  # can divide the total value by 1 
  sumA <- sum(deltaA[])
  sumB <- sum(deltaB[])
  
  #Employing "delta trick" to give vector of weights dirichlet priors
  #this is doing the dirichlet in two steps 
  #see Ogle et al. 2015 SAM model paper in Ecology Letters
  for(t in 1:n.lag){ #for total number of lags
   #the weigths for tmax - getting them to sum to 1
    wA[t] <- deltaA[t]/sumA
    #the weights for ppt - summing to 1
    wB[t] <- deltaB[t]/sumB
    #both follow relatively uninformative gamma priors
    deltaA[t] ~ dgamma(1,1)
    deltaB[t] ~ dgamma(1,1)
  }
  
  # Hyperpriors (community level)
  #easier to make uniform (~beta(1,1)) prior on non-logistic scale
  a0.mean ~ dbeta(1,1)
  mu.a0 <- logit(a0.mean)
  sd.a0 ~ dunif(0, 5)
  tau.a0 <- pow(sd.a0, -2)

  mu.tant ~ dunif(-5,5)
  sd.tant ~ dunif(0,5)
  tau.tant <- pow(sd.tant,-2)

  mu.ppt ~ dunif(-5,5)
  sd.ppt ~ dunif(0,5)
  tau.ppt <- pow(sd.ppt,-2)

  }
  