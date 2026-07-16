MoransIpermutation = function(X, Y, Cx, Ey, etas = eval(formals(sbivarSingle)$etas),
                              nPerms, normX = "none", normY = "none", pseudoCount = 1e-8, verbose = TRUE){
    if(verbose){
        cat("Starting Moran's I permutation")
    }
    if(!exists("moveTwoCoords")){
        moveTwoCoords = sbivar:::moveTwoCoords
        buildWeightMat = sbivar:::buildWeightMat
        sbivarSingle = sbivar:::sbivarSingle
        CCT = sbivar:::CCT
    }
    rownames(Cx) = rownames(X);rownames(Ey) = rownames(Y)
    movedCoords = moveTwoCoords(Cx, Ey)
    Cx = movedCoords$Cx;Ey = movedCoords$Ey
    X = normMat(X, normX, pseudoCount);Cx = Cx[rownames(X),]
    Y = normMat(Y, normY, pseudoCount);Ey = Ey[rownames(Y),]
    #No need for matching, also works disjointly
    allPvals = vapply(etas, FUN.VALUE = matrix(0, ncol(X), ncol(Y)), function(e){
        wm = buildWeightMat(Cx, Ey, wo = "Gauss", eta = e)
        Ixy = crossprod(X, wm %*% Y)
        permIxy = simplify2array(loadBalanceBplapply(integer(nPerms), function(j){
            crossprod(X[sample(nrow(X)), ], wm %*% Y[sample(nrow(Y)),])
        }))
        getPermPvalMat(Ixy, permIxy)
    })
    apply(allPvals, c(1,2), CCT)
}
CCtestEtas = function(X, Y, Cx, Ey, nPerms, wo = "Gauss", etas = eval(formals(sbivarSingle)$etas),
                      normX = "none", normY = "none", pseudoCount = 1e-8, verbose = TRUE){
    if(verbose){
        cat("Starting Moran's I random shift")
    }
    if(!exists("moveTwoCoords")){
        moveTwoCoords = sbivar:::moveTwoCoords
        buildWeightMat = sbivar:::buildWeightMat
        sbivarSingle = sbivar:::sbivarSingle
        CCT = sbivar:::CCT
    }
    rownames(Cx) = rownames(X);rownames(Ey) = rownames(Y)
    movedCoords = moveTwoCoords(Cx, Ey)
    Cx = movedCoords$Cx;Ey = movedCoords$Ey
    X = normMat(X, normX, pseudoCount);Cx = Cx[rownames(X),]
    Y = normMat(Y, normY, pseudoCount);Ey = Ey[rownames(Y),]
    allPvals = vapply(etas, FUN.VALUE = matrix(0, ncol(X), ncol(Y)), function(e){
        wm = buildWeightMat(Cx, Ey, wo = wo, eta = e)
        CC.test.Moran(X = X, Y = Y, Cx = Cx, Ey = Ey, wm = wm,
                      N.shifts = nPerms, wo = wo, eta = e)
    })
    if(dim(allPvals)[3]==1)
        allPvals[,,1]
    else
        apply(allPvals, c(1,2), CCT)
}
