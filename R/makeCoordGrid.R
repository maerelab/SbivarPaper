makeCoordGrid = function(n, Min, Max){
    xSeq = seq(Min, Max, length.out = sqrt(n))
    mat = cbind("x" = rep(xSeq, length.out = n),
                "y" = rep(xSeq, each = sqrt(n)))
    rownames(mat) = seq_len(n)
    mat
}