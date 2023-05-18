model{
  for(s in 1:n.species){
    for(i in 1:n.sites){
      y[s,i] ~  dbern(p[s,i])
      
      
      logit(p[s,i]) <- b0[s] + bTmax[GroupID[s]]*TmaxAnt[s,i]
      
      
      #summing the antecedent values
      TmaxAnt[GroupID[s],i] <- sum(TmaxTemp[GroupID[s],i,]) #summing across the total number of antecedent months
      #Generating each month's weight to sum above
      for(t in 1:n.lag){ #number of time steps we're going back in the past
        TmaxTemp[GroupID[s],i,t] <- Tmax[i,t]*wA[GroupID[s],t] 
    }
    }
  }

    # ANTECEDENT CLIMATE PRIORS
    #Sum of the weights for climate lag
    sumA <- sum(deltaA[]) #all the Tmax weights
    #Employing "delta trick" to give vector of weights dirichlet priors
    #this is doing the dirichlet in two steps 
    #see Ogle et al. 2015 SAM model paper in Ecology Letters
    for(s in 1:n.species){
    for(t in 1:n.lag){ #for the total number of lags
      #the weights for tmax - getting the weights to sum to 1
      wA[GroupID[s],t] <- deltaA[GroupID[s],t]/sumA
      #and follow a relatively uninformative gamma prior
      deltaA[GroupID[s],t] ~ dgamma(1,1)
    }
      }
    
  
}