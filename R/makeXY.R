# Create x and y variables with spatial structure
makeXY = function(xMat, yMat, nSims, distMat = as.matrix(stats::dist(rbind(xMat, yMat))),
                  what, boundsBetas = c(0.75, 1), numZones = 6, nuggets,
                  ranges, sigmas = rep(1,3), facNN = 0.25, widthStreak = 0.15,
    betaZonesX = NULL, betaZonesY = NULL, betas = NULL,
     scale = FALSE, streakFac = 1, scaleZone = .05,
    outcomeDistrs = rep("gaussian", 2),
    NBdispVec, scaleVec, nbLibs, gammaLibs){
    what = match.arg(what, choices = spatPats)
    n = nrow(xMat);m = nrow(yMat)
    idX = seq_len(n);idY = seq_len(m)+n
    if(is.null(betas))
        betas = matrix(runif(2*nSims, boundsBetas[1], boundsBetas[2]), nSims, 2)
    if(!missing(NBdispVec)){
        NBdispVec = sample(NBdispVec, nSims)
    }
    if(!missing(scaleVec)){
        scaleVec = sample(scaleVec, nSims)
    }
    if(!all(outcomeDistrs=="gaussian")){
        libSize1 = nbLibs
        libSize2 = gammaLibs
    }
    if(grepl("GPs", what)){
        if(what=="dependentGPsnoSAC"){
            covMatX = diag(n);covMatY = diag(m)
        } else {
            covMatX = sigmas[1]*(diag(nuggets[1], n) + expKernelNugget(distMat[idX, idX], ranges[1], nuggets[1]))
            covMatY = sigmas[2]*(diag(nuggets[2], m) + expKernelNugget(distMat[idY, idY], ranges[2], nuggets[2]))
        }
    } else {covMatX=covMatY=NULL}
    if(what=="none"){
        x = makeX(xMat, nSims, what = "independent", outcomeDistr = outcomeDistrs[1],
                  NBdispVec = NBdispVec, scaleVec = scaleVec, libSizes = libSize1)
        y = makeX(yMat, nSims, what = "independent", outcomeDistr = outcomeDistrs[2],
                  NBdispVec = NBdispVec, scaleVec = scaleVec, libSizes = libSize2)
    } else if(grepl("Gradient", what)){
        x = makeX(xMat, nSims, what = "gradient", coordId = "x", beta = betas[,1], libSizes = libSize1,
                  outcomeDistr = outcomeDistrs[1], NBdispVec = NBdispVec, scaleVec = scaleVec)
        y = switch(what,
            "orthogonalGradient" = makeX(yMat, libSizes = libSize2, nSims, what = "gradient", outcomeDistr = outcomeDistrs[2],
                    coordId = "y", beta = betas[,2], NBdispVec = NBdispVec, scaleVec = scaleVec),
            "sharedGradient" = makeX(yMat, libSizes = libSize2, nSims, what = "gradient", outcomeDistr = outcomeDistrs[2], coordId = "x",
                    beta = betas[,2], NBdispVec = NBdispVec, scaleVec = scaleVec),
            "oneGradient" = makeX(yMat, libSizes = libSize2, nSims, outcomeDistr = outcomeDistrs[2], what = "independent",
                    NBdispVec = NBdispVec, scaleVec = scaleVec),
            "oppositeGradient" = makeX(yMat, libSizes = libSize2, nSims, outcomeDistr = outcomeDistrs[2], what = "gradient",
                    coordId = "x", beta = -betas[,2], NBdispVec = NBdispVec, scaleVec = scaleVec))
    } else if(what %in% c("independentGPs", "dependentGPs", "dependentGPsNeg", "dependentGPsnoSAC")){
        if(what=="independentGPs"){
            x = makeX(xMat, nSims, what = "GP", libSizes = libSize1, covMat = covMatX, outcomeDistr = outcomeDistrs[1], NBdispVec = NBdispVec, scaleVec = scaleVec)
            y = makeX(yMat, nSims, what = "GP", libSizes = libSize2, covMat = covMatY, outcomeDistr = outcomeDistrs[2], NBdispVec = NBdispVec, scaleVec = scaleVec)
        } else if(what %in% c("dependentGPs", "dependentGPsNeg", "dependentGPsnoSAC")){
            sigMat = as.matrix(bdiag(covMatX, covMatY))
            sigMat[idX, idY] = expKernelNugget(distMat[idX, idY], ranges[3], nugget = 0)*sigmas[3]
        if(what=="dependentGPsNeg"){
                sigMat[idX, idY] = -sigMat[idX, idY]
            }
            sigMat[idY, idX] = t(sigMat[idX, idY]) #Symmetrize
            diag(sigMat) = diag(sigMat) + sigmas[3] #Add random effect variance
            normMat = t(mvtnorm::rmvnorm(nSims, mean = rep(0, n+m), sigma = sigMat))
            outList = mapply(list("x" = idX, "y" = idY), outcomeDistrs, list("x" = xMat, "y" = yMat), SIMPLIFY = FALSE,
                             FUN = function(id, od, co) {
                makeX(normMat = normMat[id, ], xMat = co, nSims = nSims, what = "GP", libSizes = if(od==outcomeDistrs[1]) libSize1 else libSize2,
                      outcomeDistr = od, NBdispVec = NBdispVec, scaleVec = scaleVec)
            })
            x = outList[[1]];y = outList[[2]]
        }
    } else if(grepl("Streak", what)){
        x = makeX(xMat, nSims, what = "streak", coordId = "x", beta = betas[, 1]*streakFac, widthStreak = widthStreak,
                  libSizes = libSize1, covMat = covMatX, outcomeDistr = outcomeDistrs[1], NBdispVec = NBdispVec, scaleVec = scaleVec)
        y = if(grepl("independentStreak", what)){
            makeX(yMat, libSizes = libSize2, nSims, covMat = covMatY, what = "streak",widthStreak = widthStreak,
                  coordId = "y", beta = betas[, 2]*streakFac, outcomeDistr = outcomeDistrs[2], NBdispVec = NBdispVec, scaleVec = scaleVec)
        } else if(grepl("dependentStreak", what)){
            makeX(yMat, libSizes = libSize2, nSims, covMat = covMatY, what = "streak", coordId = "x",
                  beta = betas[, 2]*streakFac, outcomeDistr = outcomeDistrs[2], NBdispVec = NBdispVec,
                  scaleVec = scaleVec, widthStreak = widthStreak)
        } else if(what == "oppositeStreak"){
            makeX(yMat, libSizes = libSize2, nSims, what = "streak",covMat = covMatY, coordId = "x",
                  beta = -betas[, 2]*streakFac, outcomeDistr = outcomeDistrs[2], NBdispVec = NBdispVec,
                  scaleVec = scaleVec, widthStreak = widthStreak)
        }
    } else if(grepl("Zones", what)){
        #Take window for signal wider
        Win = owin(range(xMat[, "x"]) + c(-0.05, 0.05),range(xMat[, "y"]) + c(-0.05, 0.05))
        xLocations = coords(runifpoint(numZones, Win))
        #xLocations = expand.grid("x" = seq(0.2, 0.8, length.out = sqrt(numZones)), "y" = seq(0.2, 0.8, length.out = sqrt(numZones)))
        if(is.null(betaZonesX))
           betaZonesX = matrix(runif(boundsBetas[1], boundsBetas[2], n = numZones*nSims), nSims, numZones, byrow = TRUE)
        if(is.null(betaZonesY))
           if(grepl("dependentZones", what)){
                betaZonesY = betaZonesX
                yLocations = xLocations
           } else if(grepl("independentZones", what)){
                betaZonesY = betaZonesX
                id <- sample(numZones, numZones/2);betaZonesY[,id] = -betaZonesY[,id]
                yLocations = xLocations
           }
        x = makeX(xMat, nSims, what = "zones", covMat = covMatX, libSizes = libSize1, beta = betaZonesX, zoneLocations = xLocations,
                  outcomeDistr = outcomeDistrs[1], NBdispVec = NBdispVec, scaleVec = scaleVec, scaleZone = scaleZone)
        y = makeX(yMat, nSims, what = "zones", covMat = covMatY, outcomeDistr = outcomeDistrs[2], NBdispVec = NBdispVec, scaleVec = scaleVec,
                  beta = switch(what, "oppositeZones"=-betaZonesY, betaZonesY), zoneLocations = yLocations, scaleZone = scaleZone)
    }
    if(scale){
        x = scale(x);y=scale(y)
    }
    rownames(x) = seq_len(nrow(x));rownames(y) = seq_len(nrow(y))
    colnames(x) = paste0(colnames(x), "X");colnames(y) = paste0(colnames(y), "Y")
    return(list("x" = x, "y" = y))
}
makeX = function(xMat, nSims, what = c("independent", "gradient", "GP", "streak", "zones", "concentric"),
                 beta, distMat, coordId, covMat=NULL, outcomeDistr = c("gaussian", "negbin", "gamma"), zoneLocations,
                 widthStreak, sdNormalVec = rep(1, nSims), NBdispVec = runif(nSims, .01, .2), normMat, scaleZone,
                 NBrelVec = runif(nSims, 0,1), libSizes = rpois(nrow(xMat), 1e5), scaleVec = runif(nSims, 0.25, 1)){
    what = match.arg(what);outcomeDistr = match.arg(outcomeDistr)
    n = nrow(xMat)
    mm = makeMeanMat(xMat, nSims, what = if(grepl("GPs", what)) "GP" else what, beta,
                     distMat, coordId, widthStreak, zoneLocations, scaleZone = scaleZone)
    out = if(!is.null(covMat) || !missing(normMat)){
        if(missing(normMat))
            normMat = t(mvtnorm::rmvnorm(nSims, mean = rep(0, n), sigma = covMat)) + mm
        if(outcomeDistr=="gaussian"){
            t(t(normMat)*sdNormalVec)
        } else {
            copula = pnorm(normMat)
            if(outcomeDistr=="negbin"){
                qnbinom(copula, size = rep(NBdispVec, each = n),
                        mu = makeMeanMatCount(mm, NBrelVec, libSizes))
            } else if(outcomeDistr=="gamma"){
                qgamma(copula, shape = rep(1/scaleVec, each = n),
                       scale = makeMeanMatCount(mm, NBrelVec, libSizes)*rep(scaleVec, each = n))
            }
        }
    } else {
        tmp = switch(outcomeDistr,
           "gaussian" = rnorm(n*nSims, sd = sdNormalVec, mean  = mm),
           "negbin" = rnbinom(n*nSims, size = rep(NBdispVec, each = n),
                              mu = makeMeanMatCount(mm, NBrelVec, libSizes)),
           "gamma" = rgamma(n*nSims, shape = rep(1/scaleVec, each = n),
                            scale = makeMeanMatCount(mm, NBrelVec, libSizes)*rep(scaleVec, each = n))
        )
        matrix(tmp, n, nSims, byrow = FALSE)
    }
    colnames(out) = paste0("feat", seq_len(nSims))
    rownames(out) = rownames(xMat)
    return(out)
}
makeMeanMat = function(xMat, nSims, what = c("independent", "gradient", "GP", "streak", "zones"),
                 beta, distMat, coordId, widthStreak, zoneLocations, scaleZone){
    what = match.arg(what)
    n = nrow(xMat)
    out = if(grepl("independent", what) || grepl("GP", what)){
        matrix(0, n, nSims)
    } else if(what == "gradient"){
        outer(xMat[,coordId], beta)
    } else if(what == "streak"){
        outer(exp(-((xMat[,coordId]-0.5)/widthStreak)^2), beta)
    } else if (what=="zones"){
        cDist = crossdist(xMat[, "x"], xMat[, "y"], zoneLocations[, "x"], zoneLocations[, "y"])
        tcrossprod(exp(-cDist^2/scaleZone), beta)
    }
    return(out)
}
makeMeanMatCount = function(mm, relVec, libSizes){
    tmp = t(t(exp(mm))*relVec)
    tmp/rowSums(tmp)*libSizes
}