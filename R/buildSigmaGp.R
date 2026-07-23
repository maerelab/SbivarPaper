#Build the sigma matrix for a Gaussian process
buildSigmaGp = function(pars, distMat, sparse = FALSE){
    pars["sigma"]^2*(diag(nrow(distMat))*pars["nugget"] +
                         expKernelNugget(distMat, pars["range"], nugget = pars["nugget"], sparse = sparse))
}
