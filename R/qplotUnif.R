#Build qq-plot with respect to standard uniform distribution
qplotUnif = function(pp, ...){
    plot(sort(pp[!is.na(pp)]), seq(0,1, length.out = sum(!is.na(pp))), ..., cex = 0.5,
         ylab = "Expected quantiles", xlab ="Observed p-values", xlim = c(0,1), ylim = c(0,1))
    abline(0,1, lty = "dashed")
}