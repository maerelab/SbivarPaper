corTests = function(X, Y, normX = "none", normY = "none", pseudoCount = 1e-8){
    X = normMat(X, normX, pseudoCount)
    Y = normMat(Y, normY, pseudoCount)
    sharedNames = intersect(rownames(X), rownames(Y))
    X = X[sharedNames,];Y = Y[sharedNames,]
    featGrid = expand.grid("featX" = colnames(X), "featY" = colnames(Y))
    out = simplify2array(loadBalanceBplapply(seq_len(nrow(featGrid)), function(i){
        unlist(cor.test(X[, featGrid[i, "featX"]], Y[, featGrid[i, "featY"]])[c("estimate", "p.value")])
    }))
    colnames(out) = apply(featGrid, 1, paste, collapse = "__")
    rownames(out) = c("Correlation", "pVal")
    t(out)
}