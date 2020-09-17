##################################################################
#
# This script contains:
# * an R function which finds optimal lambda,
# * an example of usage.
#
##################################################################


# function ----

freqs_shrink <- function(p, n, 
                         procedure=c("unif", "unif.se", "indep", "indep.se")){
  
  if(!(procedure %in% c("unif", "unif.se", "indep", "indep.se"))){
    stop("Procedue name must be one of: unif, unif.se, indep, indep.se.")
  }
  
  nx <- nrow(p)
  ny <- ncol(p)
  
  lambda <- NA
  p.shrink <- base::matrix(NA, nx, ny)
  
  if(procedure == "unif"){
    
    p.ml <- p
    p.unif  <- (rep(1, nx)/nx) %o% (rep(1, ny)/ny)
    
    mom2.ml <- sum(p.ml*((n - 1)*p.ml + 1)/n)
    var.ml  <- sum(p.ml*(1 - p.ml)/n)
    
    lambda <- var.ml/(mom2.ml + 1/(nx*ny) - 2*sum(p.ml)/(nx*ny))
    
    lambda <- ifelse(lambda > 1, 1, lambda)
    lambda <- ifelse(lambda < 0, 0, lambda)
    
    p.shrink <- lambda*p.unif + (1 - lambda)*p.ml
    
  }
  if(procedure == "unif.se"){
    
    p.ml <- p
    p.unif <- (rep(1, nx)/nx) %o% (rep(1, ny)/ny)
    
    lambda <- 0
    
    for(i in 1:nx){
      for(j in 1:ny){
        
        p.tmp <- p*n
        p.tmp[i, j] <- p.tmp[i, j] - 1
        p.tmp <- p.tmp/(n-1)
        
        p.ml.tmp <- p.tmp
        p.unif.tmp <- (rep(1, nx)/nx) %o% (rep(1, ny)/ny)
        
        lambda <- lambda +
          ((p.ml[i,j])^2 - p.ml[i,j]*p.unif[i,j] + p.ml[i,j]*p.unif.tmp[i,j] - p.ml[i,j]*p.ml.tmp[i,j])
      }
    }
    
    lambda <- lambda/sum((p.unif - p.ml)^2)
    
    lambda <- ifelse(lambda > 1, 1, lambda)
    lambda <- ifelse(lambda < 0, 0, lambda)
    
    p.shrink <- lambda*p.unif + (1 - lambda)*p.ml
  }
  if(procedure == "indep"){
    
    p.ind.x <- apply(p, 1, sum)
    p.ind.y <- apply(p, 2, sum)
    
    p.ml <- p
    
    p.ind     <- p.ind.x %o% p.ind.y
    p.ind.x.m <- p.ind.x %o% rep(1, length(p.ind.y))
    p.ind.y.m <- rep(1, length(p.ind.x)) %o% p.ind.y
    
    var.ml     <- sum(p.ml*(1 - p.ml)/n)
    
    mom2.ml    <- sum(p.ml*((n - 1)*p.ml + 1)/n)
    
    mom2.ind   <- sum(((n - 1)*(n - 2)*(n - 3)*(p.ind^2) +
                         (n - 1)*(n - 2)*p.ind*(p.ind.x.m + p.ind.y.m + 4*p.ml) +
                         (n - 1)*(2*p.ml*(p.ind.x.m + p.ind.y.m) +
                                    2*p.ml^2 + p.ind) + p.ml)/n^3)
    
    cov.ml.ind <- sum(p.ml*((n - 1)*(p.ind.x.m + p.ind.y.m - 2*p.ind) + 1 - p.ml)/n^2)
    
    mom.ml.ind <- sum(p.ml*((n - 1)*((n - 2)*p.ind + p.ind.x.m + p.ind.y.m + p.ml) + 1)/n^2)
    
    lambda <- (var.ml - cov.ml.ind) / (mom2.ml + mom2.ind - 2*mom.ml.ind)
    
    lambda <- ifelse(lambda > 1, 1, lambda)
    lambda <- ifelse(lambda < 0, 0, lambda)
    
    p.shrink <- lambda*p.ind + (1 - lambda)*p.ml
  }
  if(procedure == "indep.se"){
    
    p.ind.x <- apply(p, 1, sum)
    p.ind.y <- apply(p, 2, sum)
    
    p.ml <- p
    
    p.ind <- p.ind.x %o% p.ind.y
    
    lambda <- 0
    
    for(i in 1:nrow(p)){
      for(j in 1:ncol(p)){
        
        p.tmp <- p*n
        p.tmp[i, j] <- p.tmp[i, j] - 1
        p.tmp <- p.tmp/(n-1)
        
        p.ml.tmp <- p.tmp
        px.tmp <- apply(p.ml.tmp, 1, sum)
        py.tmp <- apply(p.ml.tmp, 2, sum)
        p.ind.tmp <- px.tmp %o% py.tmp
        
        lambda <- lambda +
          ((p.ml[i,j])^2 - p.ml[i,j]*p.ind[i,j] + p.ml[i,j]*p.ind.tmp[i,j] - p.ml[i,j]*p.ml.tmp[i,j])
      }
    }
    
    lambda <- lambda/sum((p.ind - p.ml)^2)
    
    lambda <- ifelse(lambda > 1, 1, lambda)
    lambda <- ifelse(lambda < 0, 0, lambda)
    
    p.shrink <- lambda*p.ind + (1 - lambda)*p.ml
  }

  return(list(p.shrink = p.shrink,
              lambda = lambda))
}


# example ----

SS <- function(p.est, p){ # sum of squares
  return(sum((p.est - p)^2))
}

set.seed(123)

n <- 100 ## observation number
k <- 10  ## number of levels of X and Y

p <- matrix(runif(100), k, k) 
p <- p/sum(p) ## matrix of probabilities: p1[i,j] = P(X=i, Y=j)

data <- expand.grid(x=1:k, y=1:k)[sample(1:n, prob=as.vector(p), replace=TRUE, size=n), ] ## generated data
p.est <- as.data.frame.matrix(table(data))/n ## ML estimator of p

res_unif <- freqs_shrink(p=p.est, n=n, procedure="unif")
res_unif$p.shrink 
res_unif$lambda
SS(res_unif$p.shrink, p)

res_unif.se <- freqs_shrink(p=p.est, n=n, procedure="unif.se")
SS(res_unif.se$p.shrink, p)

res_indep <- freqs_shrink(p=p.est, n=n, procedure="indep")
SS(res_indep$p.shrink, p)

res_indep.se <- freqs_shrink(p=p.est, n=n, procedure="indep.se")
SS(res_indep.se$p.shrink, p)

SS(p.est, p) 

results <- data.frame(lambda = c(0, 
                                 res_unif$lambda, res_unif.se$lambda, 
                                 res_indep$lambda, res_indep.se$lambda),
                      SS = c(SS(p.est, p), 
                             SS(res_unif$p.shrink, p), SS(res_unif.se$p.shrink, p), 
                             SS(res_indep$p.shrink, p), SS(res_indep.se$p.shrink, p)))
rownames(results) <- c("ML", "unif", "unif.se", "indep", "indep.se")
t(results)

## comments
# observe: SS(unif.se) < SS(unif)
#          SS(indep.se) < SS(indep)
#          SS(any procedure) < SS(ML)

