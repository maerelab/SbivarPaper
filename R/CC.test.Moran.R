CC.test.Moran = function (X, Y, Cx, Ey, wm, N.shifts = 999, radius = 0.5, ...)
{
    if(!exists("moveTwoCoords")){
        moveTwoCoords = sbivar:::moveTwoCoords
        buildWeightMat = sbivar:::buildWeightMat
        sbivarSingle = sbivar:::sbivarSingle
    }
    n = nrow(X);m = nrow(Y)
    obsIxy <- crossprod(X, wm %*% Y)
    jumps <- runifdisc(N.shifts, radius = radius)
    rangeFac = c(-1,1)*1e-6
    CxP = ppp(Cx[, "x"], Cx[, "y"], xrange = range(Cx[, "x"]) + rangeFac,
              yrange = range(Cx[, "y"]) + rangeFac)
    EyP = ppp(Ey[, "x"], Ey[, "y"], xrange = range(Ey[, "x"]) + rangeFac,
              yrange = range(Ey[, "y"]) + rangeFac)
    unWin = union.owin(CxP$window,EyP$window)
    simIxy = lapply(seq_len(N.shifts), function(ii){#loadBalanceBp
        # Correction "variance"
        jump = coords(jumps)[ii,]
        #Find reduced window Wc
        W.reduced <- intersect.owin(unWin, shift(CxP$window, jump))
        #Subset X in Wc, and shift back
        idinX = inside.owin(CxP, w = W.reduced)
        CxP.reduced <- shift(CxP[idinX], -jump)
        #Subset Y in Wc shifted back, and keep its coordinates
        idinY = inside.owin(EyP, w = shift(W.reduced, -jump))
        EyP.reduced <- EyP[idinY]
        wm = buildWeightMat(coords(CxP.reduced), coords(EyP.reduced), ...)
        Ixysim <- crossprod(X[idinX, ], wm %*% Y[idinY, ])
        nSim = sqrt(npoints(CxP.reduced)* npoints(EyP.reduced))
        list("Ixy" = Ixysim, "nSim" = nSim)
    })
    allStats = abind(obsIxy, vapply(simIxy, FUN.VALUE = obsIxy, function(x) x$Ixy))
    meanStats = rowMeans(allStats, dim = 2, na.rm = TRUE)
    allSamSizes = sqrt(c(sqrt(n*m), vapply(simIxy, FUN.VALUE = double(1), function(x) x$nSim)))
    allStatsNorm = (allStats-c(meanStats))*rep(allSamSizes, each = ncol(X)*ncol(Y))
    getPermPvalMat(allStatsNorm[,,1], allStatsNorm[,,-1])
}
