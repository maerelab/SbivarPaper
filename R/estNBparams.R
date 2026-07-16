estNBparams = function(x){
    Libs = rowSums(x)
    x = x[id <- (Libs > 0),]
    apply(x[, colSums(x)>0], 2, function(y){
        nbFit = glm.nb2(y = y, reg = matrix(1, nrow(x),1), s = Libs[id])
        c(theta = nbFit$theta, coef = nbFit$betas)
    })
}
Rmvnegbin = function (n, mu, Sigma, ks, empirical, ...)
{
    Cor <- cov2cor(Sigma)
    if (missing(mu))
        stop("mu is required")
    if (dim(mu)[2] != dim(Sigma)[2])
        stop("Sigma and mu dimensions don't match")
    d <- dim(mu)[2]
    normd <- MASS::mvrnorm(n, mu = rep(0, d), Sigma = Cor, empirical = empirical) #The normal-to-anything framework
    unif <- pnorm(normd)
    data <- t(qnbinom(t(unif), mu = t(mu), size = ks, ...))
    data <- .fixInf(data)
    return(data)
}

##An auxiliary function
.fixInf <- function(data) {
    # hacky way of replacing infinite values with the col max + 1
    if (any(is.infinite(data))) {
        data <-  apply(data, 2, function(x) {
            if (any(is.infinite(x))) {
                x[ind<-which(is.infinite(x))] <- NA
                x[ind] <- max(x, na.rm=TRUE)+1
            }
            x
        })
    }
    data
}
glm.nb2 = function(y, reg, s, maxit = 200L, convTol = 1e-4,
                   betas = c(log(mean(y/s)), rep(1e-10, ncol(reg)-1)),
                   theta = if(single) 0.1 else rep(0.1, length(y)),
                   single = TRUE, xFac = NULL, se = TRUE, ...){
    iter = 1L
    convergence = FALSE
    foo = try(silent = TRUE,
              while(!convergence & (iter<= maxit)){
                  betasOld = betas; thetasOld = theta
                  betas = nleqslv(x = betas, fn = ScoreNB, jac = JacNB, reg = reg , y = y, thetas = theta, s = s, ...)$x
                  mu = exp(reg %*% betas)*s
                  if(iter==1) {
                      theta =  if(single) {
                          max(1e-7,theta.mm2(y, mu))
                      } else {
                          tapply(seq_along(y), xFac, function(i){theta.mm2(y = y[i], mu = mu[i])})
                      }
                  } #Start with MoM estimate
                  theta = if(single) {
                      theta.ml2(theta = theta, y = y, mu = mu)
                  } else {
                      tapply(seq_along(y), xFac,
                             function(i){theta.ml2(theta = theta[i][1],y = y[i], mu = mu[i])})[xFac]
                  }
                  iter = iter + 1L
                  convergence = (all(abs(betas-betasOld) < convTol)) & all(abs(theta-thetasOld) < convTol)
              })
    if(inherits(foo,"try-error")){
        return(list(betas = rep(NA, ncol(reg)), theta =NA, vcov = matrix(NA, ncol(reg), ncol(reg))))
    }
    if(!convergence){warning("No convergence achieved!\n")}
    return(list(betas = betas, theta = theta))
}
theta.ml2 = function (y, mu, theta,...)
{
    nleqslv(theta, fn = scoreOD, jac = infoOD, y = y, mu = mu, ...)$x
}
scoreOD <- function(th, mu, y) {
    sum( (digamma(y+th) -digamma(th) + log(th) + 1 - log(th + mu) - (y +th)/(mu + th)))
}
infoOD <- function(th, mu, y, ...){
    -sum((-trigamma(th +y) + trigamma(th) - 1/th + 2/(mu + th) - (y + th)/(mu +th)^2))
}
theta.mm2 = function(y, mu){
    length(y)/sum((y/mu - 1)^2)
}
ScoreNB = function(betas, reg, y, thetas, s){
    mu = c(exp(reg %*% betas)*s)
    crossprod(reg,((y-mu)/(1+mu/thetas)))
}
JacNB = function(betas, reg, y, thetas, s){
    mu = c(exp(reg %*% betas)*s)
    -crossprod(reg*c((1+y/thetas)*mu/(1+mu/thetas)^2), reg)
}