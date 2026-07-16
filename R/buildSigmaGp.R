#Build the sigma matrix for a Gaussian process
buildSigmaGp = function(pars, distMat, sparse = FALSE){
    pars["sigma"]^2*(diag(nrow(distMat))*pars["nugget"] +
                         expKernelNugget(distMat, pars["range"], nugget = pars["nugget"], sparse = sparse))
}
#Optimize the gp
objGp = function(par, x, y, regMat, distMat, full, noY){
    meanVec = regMat %*% par[seq_len(NCOL(regMat))]
    sigMat = buildSigmaGp(par, distMat, n = length(x), m = length(y), noY = noY,
                          full = full, numBetas = NCOL(regMat))
    -mvtnorm::dmvnorm(if(noY) x else c(x,y), meanVec, sigMat, log = TRUE)
}
optimize_gp <- function(x, y, distMat, full = TRUE, maxEval = 1e3, init_params, regMat) {
    md <- min(distMat[distMat!=0])
    noY = missing(y)
    if(missing(init_params)){
        if(noY){#Single variable gp
            init_params = c(mean(x), sd(x), log(md*10), logit(0.1))
        } else  {
            init_params = c(mean(x), mean(x) + mean(y), sd(x), sd(y), rep(log(md/100), 2 + full))
        }
    }
    res = nloptr(x0 = init_params, eval_f = objGp,
                     lb = rep(-Inf, length(init_params)), # small, nonzero lower bound
                     ub = rep(Inf, length(init_params)), # No upper bound
                     opts = list("algorithm" = "NLOPT_LN_NELDERMEAD", "maxeval" = maxEval, "xtol_rel" = 1e-8),
                 full = full, x = x, y = if(noY) 0 else y, noY = noY, distMat = distMat, regMat = regMat
        )
    return(res)
}
#Old LRT code
# sol = optimize_gp(x, y, distMat, regMat = regMat,  init_params = c(solF$solution, -10))
# solF = optimize_gp(x, y, distMat, regMat = regMat, full = FALSE,
#                    init_params = sol$solution[-length(sol$solution)])
# #The likelihood ratio test. Be careful, the null is on the boundary of the parameter space
# # The null distribution is an even mixture of the chi-squared statistic and a point mass at 0
# # So divide the p-value by 2
# pVal = pchisq(df = 1, chiStat <- 2*(solF$objective - sol$objective), lower.tail = FALSE)/2
# #The covariance range
# range = exp(sol$solution[length(sol$solution)])

