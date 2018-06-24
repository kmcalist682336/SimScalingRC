sample.bb <- function(aa, bb, rr, mm, sig, d, L){
  new.bb <- bb
  aarr <- (aa^2)*rr
  aarr.l <- apply(aarr,2,sum)
  aarr.l <- matrix(nrow = d, rep(aarr.l,d))
  sig.dl <- 1/(sig + aarr.l)
  for(i in 1:L){
    mdj <- mm - (new.bb[,-i]%*%t(rr[,-i]*aa[,-i]))
    aarr.i <- aa[,i]*rr[,i]
    bb.draw <- c()
    for(j in 1:d){
      mu.dl <- sig.dl[j,i]*sum(aarr.i*mdj[j,])
      bb.draw[j] <- rnorm(1,mu.dl,sqrt(sig.dl[j,i]))
    }
    new.bb[,i] <- bb.draw
  }
  nbb.norm <- matrix(ncol = L, rep(sqrt(1 + apply(new.bb^2,1,sum)),L))
  new.bb <- new.bb/nbb.norm
  return(new.bb)
}

sample.sig <- function(bb,sig,d,L){
  new.sig <- sig
  bb2 <- bb^2
  for(i in 1:d){
    for(j in 1:L){
      new.sig[i,j] <- rgamma(1,.5,.5*bb2[i,j])
    }
  }
  return(new.sig)
}

sample.rr <- function(aa,bb,rr,mm,eta.l,p,L){
  new.rr <- rr
  for(j in 1:L){
    mdl <- mm - (bb[,-j]%*%t(new.rr[,-j]*aa[,-j]))
    bpb <- t(as.matrix(bb[,j]))%*%as.matrix(bb[,j])
    lel <- log(eta.l[j])
    l1mel <- log(1 - eta.l[j])
    for(i in 1:p){
      p1 <- lel - (.5*((bpb*aa[i,j]) - (2*t(as.matrix(bb[,j]))%*%as.matrix(mdl[,i])*aa[i,j])))
      p0 <- l1mel
      pr <- exp(p1)/(exp(p0) + exp(p1))
      pr <- max(pr,.001)
      pr <- min(pr,.999)
      new.rr[i,j] <- rbinom(1,1,pr)
    }
  }
  return(new.rr)
}

sample.eta.l <- function(rr,p,L){
  new.eta.l <- eta.l
  s.rr <- apply(rr,2,sum)
  s.1r <- p - s.rr
  for(i in 1:L){
    b1 <- (1/1000) + s.rr[i]
    b2 <- (999/1000) + s.1r[i]
    new.eta.l[i] <- rbeta(1,b1,b2)
  }
  return(new.eta.l)
}
