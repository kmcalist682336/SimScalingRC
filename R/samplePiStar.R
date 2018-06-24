#Sample the local level buffet parameters, pi.star
#This is the most computationally intensive draw in the sampler
#Maybe improve later
#Draw by dimensions, multivariate slice sampling
#100 draws per multiplier, reduce by order of .1 each time there is a complete failure
#No evidence that this takes more than 100 draws most of the time
#Arguments: pi.star is the p x K matrix of current pi.star values
#Kappa is the p x p matrix of normalized distances
#zz is the p x K binary matrix of features
#eta.k is the k-vector of global level buffet parameters
#p is the number of items
#K is the current number of dimensions
sample.pi.star <- function(pi.star,kappa,zz,eta.k,p,K){
  new.pi.star <- pi.star
  for(j in 1:K){
    new.pis.j <- kappa%*%new.pi.star[,j]
    #Compute matrix of leftovers
    #pi_{row,k} - (kappa_{row,column}*pi.star_{column,k})
    pi.leftovers <- matrix(ncol = p, nrow = p)
    for(i in 1:p){
      for(k in 1:p){
        pi.leftovers[i,k] <- new.pis.j[i] - (kappa[i,k]*new.pi.star[k,j])
      }
    }
    min.bumps <- pi.leftovers
    for(i in 1:p){
      for(k in 1:p){
        min.bumps[i,k] <- pi.leftovers[i,k] + (kappa[i,k]*(1 - zz[i,j]))
      }
    }
    #Compute current likelihoods
    curr.liks <- dbinom(zz[,j],1,new.pis.j)
    curr.liks <- matrix(ncol = p, nrow = p, rep(curr.liks,p), byrow = F)
    min.liks <- matrix(ncol = p, nrow = p)
    for(i in 1:p){
      for(k in 1:p){
        min.liks[i,k] <- dbinom(zz[i,j],1,min.bumps[i,k])
      }
    }
    lik.range <- abind(min.liks,curr.liks,along = 3)
    draw.u <- function(x){
      return(runif(1,x[1],x[2]))
    }
    aux.vars <- apply(lik.range,c(1,2),draw.u)
    zz.mat <- matrix(ncol = p, nrow = p, rep(zz[,j],p), byrow = F)
    inv.range <- abind(aux.vars,zz.mat,kappa,pi.leftovers,along = 3)
    inv.u <- function(x){
      if(x[2] == 1){
        return((x[1] - x[4])/x[3])
      }else{
        return((x[1] + x[4] - 1)/(-x[3]))
      }
    }
    r1 <- apply(inv.range,c(1,2),inv.u)
    rr1 <- abind(r1,zz.mat,along = 3)
    rr.min <- apply(rr1,c(1,2),min)
    rr.max <- apply(rr1,c(1,2),max)
    range.min <- apply(rr.min,2,max)
    range.max <- apply(rr.max,2,min)
    curr.beta.ll <- sum(dbeta(new.pi.star[,j],eta.k[j],1-eta.k[j], log = T))
    aux.beta.ll <- curr.beta.ll - rexp(1,1)
    st <- 0
    cc <- 1.1
    beta.eval <- function(x,off){
      gg <- sum(dbeta(x,eta.k[j],1-eta.k[j], log = T))
      if(gg > off){
        return(1)
      }else{
        return(0)
      }
    }
    while(st == 0){
      cc <- cc - .1
      rmn <- new.pi.star[,j] - (cc*(new.pi.star[,j] - range.min))
      rmx <- new.pi.star[,j] + (cc*(range.max - new.pi.star[,j]))
      props <- matrix(ncol = p, nrow = 100)
      for(i in 1:p){
        props[,i] <- runif(100,cc*range.min[i],cc*range.max[i])
      }
      eprops <- apply(props,1,beta.eval,off = aux.beta.ll)
      if(sum(eprops != 0)){
        acc <- sample(which(eprops == 1),1)
        new.pi.star[,j] <- props[acc,]
        st <- 1
      }else{
        if(cc == .1){
          st <- 1
          print(j)
        }
      }
    }
    #print(j)
  }
  return(new.pi.star)
}
