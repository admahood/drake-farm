model{
    for(k in 1:n.species){  # Loop through species
      for(i in 1:n.sites){ #loop through the sites
      # Likelihood
        # Data are around occupancy probability
        #this is "apparent" occupancy because it's not
        # corrected for detection
        y[k,i] ~ dnorm(eta[k,i], sig[k])
        
        #occupancy dependent on covariates
        #Probit model for occupancy as opposed to logistic
        eta[k,i] <- a0[GroupID[k]] + 
          a1[GroupID[k]]*TAnt[i] + 
          a2[GroupID[k]]*PAnt[i] +
          inprod(lv.coef[k,], LV[i,])
        
        
        #could add in random effects and other covariates to above
        # if desired
        #could also make the a1, a2 or other parameters dependent on
        # a speciesID or species group ID (e.g. effect of temperature
        # on grasses vs. forbes, large vs. small seeds, etc)
        
      } #n.sites
      
      #new variance for model that allows residual variance
      #to be correct
      sig[k] <- (1/(1-sum(lv.coef[k,1:n.lv]^2)))

    } #n.species
  
  #Species-group level Priors
  for(g in 1:n.groups){
    a0[g] ~ dnorm(mu.a0, tau.a0)
    a1[g] ~ dnorm(mu.tant, tau.tant)
    a2[g] ~ dnorm(mu.ppt, tau.ppt)
  }
  
  #Antecedent climate values
  for(i in 1:n.sites){
    #antecedent climate value is the sum of all the
    #weighted lagged values
    TAnt[i] <- sum(TmaxTemp[i,])
    PAnt[i] <- sum(PPTTemp[i,])

    #these are where the data you input into the model
    #come in for temp and ppt - here they're multiplied
    #by the weights the model has given them for their
    #relative importance in driving the climate effect
    for(t in 1:n.lag){
      TmaxTemp[i,t] <- Tmax[i,t]*wA[t]
      PPTTemp[i,t] <- PPT[i,t]*wB[t]
    }

    #eventually - I would like to code the weights
    #to be different per group membership, so, for example,
    #annual forbes can respond to temperature in spring
    #whereas perennial grasses can respond to temperature
    # in fall
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
  
  
  #SETTING UP COVARIANCE MATRIX
  #covariance matrix is going to be positive on the diagonal,
  #zero on the topand between -1 and 1 correlatoin on the
  # bottom 
  #Latent variables
  #Re-read these two papers to remind myself the 
  #why of latent variable models
  #they can be thought of as site-level variation
  #https://esajournals.onlinelibrary.wiley.com/doi/full/10.1002/ecy.2754
  #https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12236
  
  for(j in 1:n.sites){
    for(v in 1:(n.lv)){ #number of latent variables is species/2
      #latent variable for site j and latent variable v
      LV[j,v] ~ dnorm(0,1)
    }
  }
  
  #Latent variable coefficients with constranits
  #diagonal elements positive, upper diagnoal zeros
  for(v in 1:(n.lv-1)){
    for(v2 in (v+1):n.lv){
      #set the upper diagonal to 0
      lv.coef[v, v2] <- 0
    }
  }
  
  #sign contraints on diagonal elements
  for(v in 1:n.lv){
    lv.coef[v,v] ~ dunif(0,1)
  }
  
  #lower diagonal can be any value between total
  #negative and total positive correlation
  for(v in 2:n.lv){
    for(v2 in 1:(v-1)){
      lv.coef[v, v2] ~ dunif(-1,1)
    }
  }
  
  #other elements free
  for(i in (n.lv+1):n.species){
    for(v in 1:n.lv){
      lv.coef[i,v] ~ dunif(-1,1) 
    }
  }

  }
  