plotSplines = function(x, y, xMat, yMat, newGrid, addTitle = TRUE,  k = -1, n_points = 5e2,
                       fam1 = gaussian(), fam2 = gaussian(), offset1 = NULL, offset2 = NULL,
                       scaleFun = scaleMinusOne){
    if(missing(newGrid)){
        newGrid = sbivar:::buildNewGrid(xMat, yMat, n_points = n_points)
    }
    modelx <- fitGAM(df = data.frame("value" = x, xMat), family = fam1, offset = offset1)
    modely <- fitGAM(df = data.frame("value" = y, yMat), family = fam2, offset = offset2)
    predx = vcovPredGam(modelx, newdata = newGrid)
    predy = vcovPredGam(modely, newdata = newGrid)
    corContr = (cen1 <- (predx$pred-mean(predx$pred)))*(cen2 <- (predy$pred-mean(predy$pred)))
    covEst = sum(corContr)/((nrow(newGrid)-1)*sd(predx$pred)*sd(predy$pred))
    dat = rbind(data.frame(newGrid, value = scaleFun(predx$pred), feature = "x"),
                data.frame(newGrid, value = scaleFun(predy$pred), feature = "y"),
                data.frame(newGrid, value = scaleFun(corContr), feature = "cor"))
    gridMolt = melt(dat, id.vars = c("x", "y", "feature"), value.name = "Value")
    gridMolt$feature = factor(gridMolt$feature, levels =c("x", "y", "cor"),
                              labels = c("x", "y", "cor"), ordered = TRUE)
    ggplot(gridMolt, aes(x, y, fill = Value)) +
        geom_raster() + coord_fixed() +
        facet_grid( ~ feature) +
        scale_fill_viridis_c(option = "H", name = "") +
        if(addTitle)
            labs(title = paste("Estimated spline surfaces, and contributions to correlation estimate",round(covEst, 3)))
}
plotSplinesVicari = function(gene, met, sample, ...){
    plotSplines(exprObjs[[sample]][, gene], msiObjs[[sample]][, met],
                magpieCoords[[sample]]$Stx, magpieCoords[[sample]]$Msi, newGrids[[sample]],
                fam1 = nb(), fam2 = stats::Gamma(link = "log"),
                offset1 = rowSums(exprObjs[[sample]]), offset2 = rowSums(msiObjs[[sample]]), ...)
}