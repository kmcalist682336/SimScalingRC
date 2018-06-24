#A function for approximating the log likelihood for a new proposal of A
approx.ll <- function(new.aa,zz,pi.star){
  new.dists <- make.norm.dist.mat(new.aa)
  new.pis <- new.dists%*%pi.star
  new.probs <- (zz*new.pis) + ((1 - zz)*(1 - new.pis))
  return(sum(log(new.probs)))
}
#A function for sampling AA, the matrix of latent variables in the second layer of the model
#The key innovation in this model.
#Pretty fast right now.  I spent a lot of time on this function, so I'm fine with it for now.
#Arguments: aa is p x L matrix of latent vars
#bb is d x L matrix of sedon lovel loadings
#rr is p x L binary matrix for text features
#mm is d x p matrix of augmented text data
#zz is p x K binary matrix of top level features
#pi.star is p x K matrix of top level local feature probabilities
#p is number of items
#L is the current number of features on the lower level
sample.aa <- function(aa,bb,rr,mm,zz,pi.star,p,L){
  #try an element by element slice sampler
  #set.seed(1234)
  new.aa <- aa
  for(i in 1:L){
    mdj <- mm - (bb[,-i]%*%t(rr[,-i]*new.aa[,-i]))
    b2 <- apply(bb^2,2,sum)[i]
    uc.var <- 1/(1 + (rr[,i]*b2))
    bl <- matrix(ncol = p,rep(bb[,i],p))
    mbl <- bl*mdj
    mbl <- apply(mbl,2,sum)
    uc.mean <- uc.var*mbl*rr[,i]
    curr.vals <- new.aa[,i]
    curr.ll <- dnorm(curr.vals, uc.mean, sqrt(uc.var), log = T)
    aux.uc <- curr.ll - rexp(p,1)
    range.p <- uc.mean + sqrt(-uc.var*((2*aux.uc) + log(2*pi*uc.var)))
    range.m <- uc.mean - sqrt(-uc.var*((2*aux.uc) + log(2*pi*uc.var)))
    range.l <- c()
    range.u <- c()
    for(j in 1:p){
      range.l[j] <- min(range.p[j],range.m[j])
      range.u[j] <- max(range.p[j],range.m[j])
    }
    curr.bn <- approx.ll(new.aa = new.aa, zz = zz, pi.star = pi.star)
    aux.bn <- curr.bn - rexp(1,1)
    st <- 0
    cc <- 0
    while(st == 0){
      cc <- cc + 1
      naa <- new.aa
      props <- c()
      for(j in 1:p){
        props[j] <- runif(1,range.l[j],range.u[j])
      }
      naa[,i] <- props
      eval.naa <- approx.ll(new.aa = naa, zz = zz, pi.star = pi.star)
      if(eval.naa > aux.bn){
        new.aa[,i] <- props
        st <- 1
      }else{
        if(cc == 100){
          st <- 1
        }
        diffs <- sign(props - new.aa[,i])
        for(j in 1:p){
          if(diffs[j] == -1){
            range.l[j] <- props[j]
          }else{
            range.u[j] <- props[j]
          }
        }
      }
    }
    #print(i)
  }
  return(new.aa)
}
