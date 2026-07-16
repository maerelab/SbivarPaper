estGammaParams = function(X){
    Libs = rowSums(X)
    gammaFits = apply(X[Libs>0, colSums(X)>0]+1e-8, 2, function(y){
        try(silent = TRUE, glm(y ~1, offset = log(Libs[Libs>0]), family = stats::Gamma(link = "log")))
    })
    vapply(gammaFits[sapply(gammaFits, inherits, "glm")], FUN.VALUE = double(1), function(y){
        summary(y)$dispersion
    })
}