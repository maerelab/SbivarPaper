plotVicariAll <- function(pairs, size = 0.5, shape = 18, stripTextSizeX = 6.5, axisLegendSize = 9,
                          stripTextSizeY = stripTextSizeX, samples = sampleNamesMouse,
                          facetForm = image ~ feat, xLimShift = 0,
                          legendShift = -100) {
  splitPairs <- simplify2array(strsplit(pairs, split = "__"))
  genes <- splitPairs[1, ]
  mets <- splitPairs[2, ]
  Dfs <- mapply(genes, mets, samples, FUN = function(gene, met, sample) {
    plotDf <- makePlotSingleDf(gene, met,
      sample = sample, normXMats = exprObjsNorm, normYMats = msiObjsNorm,
      metHeader = names(unMetsIdentified)[match(met, unMetsIdentified)], geneHeader = gene
    )
    plotDf$x <- plotDf$x + (1 - max(plotDf$x)) / 2
    plotDf
  }, SIMPLIFY = FALSE)
  maxY <- max(sapply(Dfs, function(x) max(x$y)))
  Plots <- lapply(seq_along(Dfs), function(i) {
    ggplot(data = Dfs[[i]], aes(x = x, y = y, col = value)) +
      geom_point(size = size, shape = shape) +
      facet_grid(facetForm) +
      scale_colour_gradient(low = "yellow", high = "blue", name = "Outcome") +
      xlab("Coordinate 1") +
      ylab("Coordinate 2") +
      coord_fixed() +
      xlim(c(0 - xLimShift, 1 + xLimShift)) +
      ylim(c(0, maxY)) +
      theme(axis.title.x = if (i == length(Dfs)) element_text(size = axisLegendSize) else element_blank())
  })
  p <- patchwork::wrap_plots(Plots, ncol = 1, guides = "collect", axes = "collect", axis_titles = "collect") &
    theme(
      legend.position = "right", strip.clip = "on",
      legend.box.spacing = unit(1, "pt"),
      legend.margin = margin(l = legendShift, b = 5, t = 0, r = 0),
      axis.title = element_text(size = axisLegendSize),
      plot.margin = margin(t = 1.25, r = 2, b = 1.25, l = 2),
      strip.text.x = element_text(
        size = stripTextSizeX,
        margin = margin(t = 2, b = 2)
      ),
      strip.text.y = element_text(size = stripTextSizeY),
      axis.text = element_blank(), axis.ticks = element_blank()
    )
  print(p)
}
plotVicariGenesMets <- function(genes, mets, size = 0.12, shape = 18, stripTextSizeX = 7,
                                metHeader = "Metabolite", geneHeader = "Gene", axisLegendSize = 7.5,
                                stripTextSizeY = 7, samples = sampleNamesMouse, norm = "rel") {
  tmp <- lapply(samples, function(sam) {
    tmpGen <- data.frame(base::do.call(rbind, lapply(genes, function(gene) {
      val <- if (gene %in% colnames(exprObjsNorm[[sam]])) exprObjsNorm[[sam]][, gene] else rep(NA, nrow(exprObjsNorm[[sam]]))
      cbind(magpieCoords[[sam]]$Stx, "value" = scaleZeroOne(val))
    })), "feat" = rep(genes, each = nrow(magpieCoords[[sam]]$Stx)))
    tmpMet <- data.frame(base::do.call(rbind, lapply(mets, function(met) {
      val <- if (met %in% colnames(msiObjsNorm[[sam]])) msiObjsNorm[[sam]][, met] else rep(NA, nrow(msiObjsNorm[[sam]]))
      cbind(magpieCoords[[sam]]$Msi, "value" = scaleZeroOne(val))
    })), "feat" = rep(names(unMetsIdentified)[match(mets, unMetsIdentified)], each = nrow(magpieCoords[[sam]]$Msi)))
    rbind(tmpGen, tmpMet)
  })
  plotDf <- data.frame(base::do.call(rbind, tmp), "image" = rep(samples, sapply(tmp, nrow)))
  plotDf$feat <- factor(plotDf$feat,
    levels = c(genes, names(unMetsIdentified)),
    labels = c(genes, names(unMetsIdentified)), ordered = TRUE
  )
  p <- ggplot(data = plotDf, aes(x = x, y = y, col = value)) +
    geom_point(size = size, shape = shape) +
    facet_grid(feat ~ image) +
    scale_colour_gradient(low = "yellow", high = "blue", name = "Outcome") +
    xlab("Coordinate 1") +
    ylab("Coordinate 2") +
    coord_fixed() +
    theme(
      legend.position = "top", axis.title = element_text(size = axisLegendSize),
      strip.text.x = element_text(size = stripTextSizeX),
      strip.text.y = element_text(size = stripTextSizeY),
      axis.text = element_blank(), axis.ticks = element_blank()
    )
  print(p)
}
plotVicariSingle <- function(gene, met, samples = sampleNamesMouse, size = 0.7,
                             shape = 18, stripTextSizeX = 8.5, axisLegendSize = 8.5, stripTextSizeY = 6,
                             facetForm = image ~ feat) {
  mh <- names(unMetsIdentified)[match(met, unMetsIdentified)]
  plotDf <- base::do.call(rbind, lapply(samples, FUN = function(sample) {
    makePlotSingleDf(gene, met,
      sample = sample, normXMats = exprObjsNorm, normYMats = msiObjsNorm,
      metHeader = mh, geneHeader = gene
    )
  }))
  plotDf$feat <- factor(plotDf$feat, levels = c(gene, mh), ordered = TRUE)
  # Add gene and metabolite in the strip titles
  p <- ggplot(data = plotDf, aes(x = x, y = y, col = value)) +
    geom_point(size = size, shape = shape) +
    facet_grid(facetForm) +
    scale_colour_gradient(low = "yellow", high = "blue", name = "Outcome") +
    xlab("Coordinate 1") +
    ylab("Coordinate 2") +
    coord_fixed() +
    theme(
      legend.position = "right", axis.title = element_text(size = axisLegendSize),
      strip.text.x = element_text(size = stripTextSizeX),
      strip.text.y = element_text(size = stripTextSizeY),
      axis.text = element_blank(), axis.ticks = element_blank()
    )
  print(p)
}
makePlotSingleDf <- function(gene, met, sample, normXMats, normYMats, metHeader = "Metabolite", geneHeader = "Gene") {
  magpieCoordsShifted <- lapply(magpieCoords, function(x) (shiftCoord(x))) # rotate_clouds
  Each <- sapply(magpieCoordsShifted[sample], function(x) nrow(x$Stx) + nrow(x$Msi))
  geneLocs <- base::do.call(rbind, lapply(magpieCoordsShifted[sample], function(x) x$Stx))
  metLocs <- base::do.call(rbind, lapply(magpieCoordsShifted[sample], function(x) x$Msi))
  geneVal <- unlist(lapply(normXMats[sample], function(xx) if (gene %in% colnames(xx)) xx[, gene] else rep(NA, nrow(xx))))
  metVal <- unlist(lapply(normYMats[sample], function(xx) if (met %in% colnames(xx)) xx[, met] else rep(NA, nrow(xx))))
  data.frame(rbind(geneLocs, metLocs),
    "image" = c(rep(sample, each = Each)),
    "feat" = factor(c(rep(geneHeader, nrow(geneLocs)), rep(metHeader, nrow(metLocs))),
      labels = c(geneHeader, metHeader), levels = c(geneHeader, metHeader), ordered = TRUE
    ),
    "value" = c(scaleZeroOne(geneVal), scaleZeroOne(metVal))
  )
}
plotGAMvicari <- function(Sam, gene, met, Families = FamiliesVicari, n_points_grid = 1e3) {
  plotGAMs(as.matrix(exprObjs[[Sam]]), msiObjs[[Sam]] + 1e-8,
    Cx = magpieCoords[[Sam]]$Stx,
    Ey = magpieCoords[[Sam]]$Msi, features = c(gene, met), families = Families,
    offsets = list(rowSums(exprObjs[[Sam]]), rowSums(msiObjs[[Sam]])), n_points_grid = n_points_grid
  )
}
plotGodfreySingle <- function(gene, met, List = godfrey, samples = names(List), pair, ...) {
  plotGodfreyAll(
    pairs = rep(if (!missing(pair)) pair else paste(gene, met, sep = "__"), length(samples)),
    samples = samples, List = List, ...
  )
}
plotGodfreySinglePair <- function(pair, samples = names(godfrey), ...) {
  plotGodfreyAll(pairs = rep(pair, length(samples)), samples = samples, ...)
}
plotGodfreyAll <- function(pairs, List = godfreyIdentified, samples = names(List),
                           facetForm = image ~ variable, stripSizeX = 5.25,
                           stripSizeY = stripSizeX, size = 0.015, axisLegendSizeX = 7.5,
                           axisLegendSizeY = axisLegendSizeX, minLibSize = 1e2, shape = 18,
                           minChar = 15, maxChar = 35, legendShift = -200, xLimShift = 0) {
  splitPairs <- simplify2array(strsplit(pairs, split = "__"))
  Max <- max(sapply(List, function(x) max(x$Coord)))
  mets <- sapply(allMetsCommon, substr, 1, maxChar)[sapply(splitPairs[2, ], function(x) grep(substr(x, 1, min(minChar, nchar(x))), names(allMetsCommon)))]
  Coords <- lapply(selfName(samples), function(Sam) {
    List[[Sam]]$Coord[, "y"] <- List[[Sam]]$Coord[, "y"] - min(List[[Sam]]$Coord[, "y"])
    List[[Sam]]$Coord <- List[[Sam]]$Coord / Max
    List[[Sam]]$Coord[, "x"] <- List[[Sam]]$Coord[, "x"] - min(List[[Sam]]$Coord[, "x"])
    List[[Sam]]$Coord
  })
  xMax <- max(sapply(Coords, function(x) max(x[, "x"])))
  yMax <- max(sapply(Coords, function(x) max(x[, "y"])))
  Coords <- lapply(Coords, function(CC) {
    CC[, "x"] <- CC[, "x"] + (xMax - max(CC[, "x"])) / 2 # Shift to center
    CC[, "y"] <- CC[, "y"] + (yMax - max(CC[, "y"])) / 2
    CC
  })
  Plots <- lapply(seq_along(samples), function(x) {
    Last <- (x == length(samples))
    Sam <- samples[x]
    dfFoo <- data.frame(
      "image" = Sam, Coords[[Sam]],
      "value" = c(
        if (splitPairs[1, x] %in% colnames(List[[Sam]]$Stx)) {
          scaleZeroOne(List[[Sam]]$Stx[, splitPairs[1, x]] / rowSums(List[[Sam]]$Stx))
        } else {
          rep(NA, nrow(List[[Sam]]$Stx))
        },
        "Metabolite" = if (length(grep(splitPairs[2, x], colnames(List[[Sam]]$Msi)))) {
          scaleZeroOne(List[[Sam]]$Msi[, grep(splitPairs[2, x], colnames(List[[Sam]]$Msi))] / rowSums(List[[Sam]]$Msi))
        } else {
          rep(NA, nrow(List[[Sam]]$Msi))
        }
      ),
      "variable" = rep(c(splitPairs[1, x], mets[x]), each = nrow(Coords[[Sam]])), row.names = NULL
    )
    ggplot(data = dfFoo, aes(x = x, y = y, col = value)) +
      geom_point(size = size, shape = shape) +
      facet_grid(facetForm) + # , space = "free", scales = "free"
      scale_colour_gradient(high = "blue", low = "yellow", name = "Outcome") +
      xlab("Coordinate 1") +
      ylab("Coordinate 2") +
      xlim(c(0 - xLimShift, xMax + xLimShift)) +
      ylim(c(0, yMax)) +
      theme(axis.title.x = if (Last) element_text(size = axisLegendSizeX) else element_blank()) +
      coord_fixed()
  })
  p <- patchwork::wrap_plots(Plots, ncol = 1, guides = "collect") &
    theme(
      legend.position = "right", legend.box.spacing = unit(1, "pt"),
      legend.margin = margin(l = legendShift, b = 5, t = 0, r = 0),
      axis.text = element_blank(), axis.ticks = element_blank(),
      axis.title.y = element_text(size = axisLegendSizeY),
      plot.margin = margin(t = 1.25, r = 1.5, b = 1.25, l = 1.5),
      strip.text.x = element_text(size = stripSizeX, margin = margin(t = 1.5, b = 1.5)),
      strip.text.y = element_text(size = stripSizeX)
    )
  print(p)
}
