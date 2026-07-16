prepareCzaplewskiWeights = function(wm){
    #Czaplewski1993
    S3 = sum(wm *t(wm))
    S4 = sum(wm^2)
    S5 = sum(wm %*% wm)
    S6 = sum(crossprod(wm, wm) + tcrossprod(wm, wm))
    S1 = S3 + S4
    S2 = 2*S5 + S6
    c("S1" = S1, "S2" = S2, "S3" = S3, "S4" = S4, "S5" = S5, "S6" = S6)
}