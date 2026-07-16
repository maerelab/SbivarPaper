#Fast geary's C implementation for permutations, as in bispdep::geary.bi
geary.bi.fast = function (varX, varY, wc, dgr)
{
    wc$n1/wc$S0 * sum(dgr * outer(varX, varY, "-")^2)#/sum(zx^2) #Scaling is pointless in permutations
}
# As in the literature
geary_bi_wartenberg <- function(x, y, L) {
    # Center variables
    x_c <- x - mean(x)
    y_c <- y - mean(y)

    # Numerator
    num <- as.numeric(crossprod(x_c, L) %*% y_c)
    #No denominator needed for permutation tests
   return(num)
}