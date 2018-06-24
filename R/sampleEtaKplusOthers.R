#Faster euclidian distance calculator
#Calculates euc distance between rows
euc_dist <- function(m){
  mtm <- Matrix::tcrossprod(m)
  sq <- rowSums(m^2)
  sqq <- outer(sq,sq,"+") - (2*mtm)
  sqq[sqq <0] <- 0
  return(sqq)
}
#Makes new matix kappa with input aa
make.norm.dist.mat <- function(m){
  t.dist <- euc_dist(abs(m))
  t.dist <- exp(-as.matrix(t.dist))
  sms <- apply(t.dist,1,sum)
  sms <- matrix(ncol = p,rep(sms,p))
  t.dist <- t.dist/sms
  return(t.dist)
}
#Makes a new p times k matrix of pis
make.new.pis <- function(kappa,pi.star){
  return(kappa%*%pi.star)
}
#Sample new values of eta.k, the global buffet parameters for the top level pi.stars
#Constructed to be sparse in pi.star
#Arguments: pi.star is p x K matrix of local prob params
#p is the number of items
#K is the current number of dimensions
sample.eta.k <- function(pi.star,p,K){
  sum.ps <- apply(pi.star,2,sum)
  new.eta.k <- c()
  for(j in 1:K){
    nek1 <- (1/1000) + sum.ps[j]
    nek2 <- (999/1000) + p - sum.ps[j]
    new.eta.k[j] <- rbeta(1,nek1,nek2)
  }
  return(as.matrix(new.eta.k))
}
