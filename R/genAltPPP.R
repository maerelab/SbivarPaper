genAltPPP = function(lambda, what, numGenes = nFeats, gpRange = 1, mu, win, scale,
                     gpRangeScale, gpRangeScaleI, p2){
    boo = lapply(seq_len(numGenes), function(i){
        ppp = switch(what,
                     "Poisson" = rpoispp(lambda),
                     "Cauchy" = rCauchy(lambda/Mu, scale = scale, mu = mu, win = win),
                     "GaussPoisson" = rGaussPoisson(lambda, gpRange*gpRangeScale, p2 = p2, win = win),
                     "MaternCluster" = rMatClust(lambda, scale = gpRange*gpRangeScale, mu = mu, win = win),
                     "MaternInhibition" = rMaternI(lambda, r = gpRange*gpRangeScaleI, win = win),
                     "Thomas" = rThomas(lambda, scale = gpRange*gpRangeScale, mu = mu, win = erosion(win, 50))
        )
        data.frame(coords(ppp), feature = paste0("gene", i))
    })
    df = do.call(what = rbind, boo)
    ppp(x = df$x, y = df$y, marks = data.frame("feature" = df$feature), drop = FALSE)
}