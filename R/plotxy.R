plotxy = function(X, Y, xlim = range(c(X[,1], Y[,1])), ylim = range(c(X[,2], Y[,2])),
    pch = 20, cex = 6, pchY = pch, xlab = "Coordinate 1", ylab = "Coordinate 2", ...){
    plot(x = as.matrix(X), xlim = xlim, ylim = ylim,
         asp = 1, pch = pch, cex = cex, xlab = xlab, ylab = ylab, ...)
    points(as.matrix(Y), col = "blue", pch = pchY, cex = cex)
}
plotx = function(X, pch = 20, cex = 6, pchY = pch, xlab = "Coordinate 1",
                 ylab = "Coordinate 2", ...){
    plot(x = as.matrix(X), xlim = range(c(X[,1])), ylim = range(c(X[,2])),
         asp = 1, pch = pch, cex = cex, xlab = xlab, ylab = ylab, ...)
}