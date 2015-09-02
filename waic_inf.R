# Function calculates waic for infection models
# input mcmc.list with posterior sims and data.frame for model input
# output list: log posterior predictive density, p_WAIC2, WAIC, mean likelihood
waic_inf <- function(posterior, data){
  with(data,{
    chains <- length(posterior)
    store <- dim(posterior[[1]])[1]
    L <- array(dim=c(N.obs, chains, store))
    L_bar <- rep(NA, N.obs)
    var_LL <- rep(NA, N.obs)
    
    for (i in 1:N.obs){
      for (j in 1:chains){
        post_sims <- posterior[[j]]
        indx_probs <- which(dimnames(post_sims)[[2]] == 
                          paste("prob[", i, "]", sep=""))
        prob_vals <- post_sims[, indx_probs]
        
        L[i, j, ] <- dbinom(rep(n.inf[i], store), 
                            size=n.recov[i], 
                            prob = prob_vals, 
                            log=TRUE)
      }
      L_bar[i] <- mean(exp(c(L[i, , ])))
      var_LL[i] <- var(c(L[i, , ]))
    }
    
    # Calculate WAIC for each virus separately
    lppd <- rep(NA, length(unique(Virus)))
    p_D2 <- rep(NA, length(unique(Virus)))
    WAIC <- rep(NA, length(unique(Virus)))
    
    for(i in 1:length(unique(Virus))){
      lppd[i] <- sum(log(L_bar[which(Virus==i)]))
      p_D2[i] <- sum(var_LL[which(Virus==i)])
      WAIC[i] <- -2 * (lppd[i] - p_D2[i])
    }
    
    return(list(lppd=lppd, p_D2=p_D2, WAIC=WAIC, L_bar))
  })
}