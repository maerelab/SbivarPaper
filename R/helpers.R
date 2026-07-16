#Some aux functions
toArray = function(x) {
    array(unlist(x[, -1]), dim = c(nrow(x), (ncol(x)-1)/2, 2))
}
#apply(y, c(1,2), var) efficiently
colVarsArray = function(x){
    rm = rowMeans(x, dims = 2)
    rowSums((x-c(rm))^2, dims = 2)/(dim(x)[3]-1)
}
#Sd without n-1
sdx = function(x) sqrt(mean((x-mean(x))^2))
#Normalize
normalize = function(x) (x-mean(x))/sdx(x)
#SumNormalize
sumnormalize = function(x) x/sum(x)
#Boxplot with diamonds for means
boxDiamond = function(x, ...){
    boxplot(x, ...)
    points(seq_len(ncol(x)), colMeans(x), pch = 5)
}
#Extract upper and lower triagonal matrices
getUpTri = function(x) x[upper.tri(x)]
getLowTri = function(x) x[lower.tri(x)]
#Get standard deviation from linear fit from lm.git
getSigma = function(Fit){
    sqrt(sum(Fit$residuals^2)/(length(Fit$residuals)-2))
}
#Predict from a linear fit
predLinFit = function(x, linFit){
   cbind(1, x)[, !is.na(linFit$coefficients)] %*% linFit$coefficients[!is.na(linFit$coefficients)]
}
#Trace function
tr = function(x, dim = c(1,2)) {
    if(is.matrix(x)|| is(x, "dgeMatrix") || is(x, "dgCMatrix")){
        sum(diag(x))
    } else if(is.array(x)){
        apply(x, dim, tr)
    } else {
        stop("Trace function not implemented for ", class(x))
    }
}
#logit
logit = function(x) log(x/(1-x))
#expit
expit = function(x) {exp(x)/(1+exp(x))}
#Histogram and qqplot
histAndQplot = function(x, main = ""){
    hist(x, main = main, xlab = "P-value")
    qplotUnif(x, main = main)
}
## ggplotSpace
plotgg = function(x, y, xMat, yMat, size = 4){
    plotDf = data.frame(rbind(xMat, yMat), "value" = c(x, y),
                        "feat" = rep(c("x", "y"), times = c(length(x), length(y))))
    ggplot(data = plotDf, aes(x=x, y=y, col = value)) + geom_point(size = size) +
        facet_grid(~feat) +
        scale_colour_gradient(low = "yellow", high = "blue", name = "Outcome") +
        xlab("Dimension 1") + ylab("Dimension 2") + coord_fixed()
}
plotggList = function(dat, xMat, yMat, i = 1, j = 1, ...){
    plotgg(dat$x[, i], dat$y[, j], xMat, yMat, ...)
}
outDiagZero = function(x){
    tmp = outer(x,x)
    diag(tmp) = 0
    tmp
}
#Convert z-value to p-value
makePval = function(z){
    z[is.na(z)] = 0
    tmp = pnorm(z, lower.tail = TRUE)
    tmp[z>0] = pnorm(z[z>0], lower.tail = FALSE)
    2*unname(tmp)
}
scaleZeroOne = function(y, na.rm = TRUE){
    (y-min(y, na.rm = na.rm))/diff(range(y, na.rm = na.rm))
}
scaleMinusOne = function(y, na.rm = TRUE){
    (y-min(y, na.rm = na.rm))/diff(range(y, na.rm = na.rm))*2-1
}
expKernel = function(distMat, range, thresh = 1e-4, sparse = FALSE){
    out = exp(-(distMat/range)^2)
    if(sparse){
        out[out<thresh] = 0
        out = as(out, "sparseMatrix")
    }
    return(out)
}
expKernelNugget = function(distMat, range, nugget, sparse = FALSE){
    (1-nugget)*expKernel(distMat, range, sparse = sparse)
}
corner = function(x, n = 6, digits = 5){
    round(x[seq_len(n), seq_len(n)], digits = digits)
}
selfName = function(x){names(x)=x;x}
logNorm = function(x, pseudoCount = 1e-8){
    log((as.matrix(x)+pseudoCount)/rowSums(x))[rowSums(x)>0,]
}
logNormScale = function(x, scale = TRUE, center = TRUE, ...){
    scale(logNorm(x, ...), scale = scale, center =center)
}
dir.nw = function(dir){
    dir.create(dir, showWarnings = FALSE)
}
getMinCombo = function(mat){
    pAdjMat = p.adjust(mat, method ="BH")
    if(all(pAdjMat > sigLevel)){
        return(NULL)
    }
    wm = which.min(mat)
    id = arrayInd(wm, dim(mat))
    met = rownames(mat)[id[1,1]]
    gene = colnames(mat)[id[1,2]]
    c('met' = met, 'gene' = gene)
}
tryMat = function(mat, met, gene){
    if(inherits(try(tmp <- mat[met, gene], silent = TRUE), "try-error")||is.null(tmp)){NA} else tmp
}
getSigPairs = function(x, collapse = TRUE){
    pAdj = p.adjust(x, method = "BH")
    id = pAdj<sigLevel
    matId = expand.grid("met" = rownames(x), "gene" = colnames(x))[id, ,drop = FALSE][ord <- order(pAdj[id]),]
    if(collapse){
        pairs = apply(matId, 1, paste, collapse = "_")
        return(pairs)
    } else return(cbind(matId, "pAdj" = pAdj[id][ord]))
}
getAllPairs = function(x, collapse = TRUE){
    mat = expand.grid("met" = rownames(x), "gene" = colnames(x))
    if(collapse){
        pairs = apply(mat, 1, paste, collapse = "_")
        return(pairs)
    } else return(mat)
}
#A is mxmxp, M is mxm
arrayMatProd = function(A, M){
    n = dim(A)[1];p = dim(A)[3]
    dn = dimnames(A)
    # # Process one slice at a time to minimize memory use
    # for (i in seq_len(p)) {
    #     A[,,i] <- M %*% A[,,i]
    # }
    # A
    # # Reshape A into n x (n*p)
    A <- matrix(A, n, n * p)
    # Multiply in one go
    result <- M %*% A
    # Reshape back to n x n x p
    result <- array(result, dim = c(n, n, p))
    dimnames(result) = dn
    return(result)
}
#A and B are (mxmxp), yields pxp
arrayProd = function(A, B){
    n = dim(A)[1];p = dim(A)[3];k = dim(B)[3]
    # Reshape each into (n^2) x p
    A_mat <- matrix(A, n*n, p)
    B_mat <- matrix(B, n*n, k)

    # p x p result: matrix of all slice inner products
    C <- crossprod(A_mat, B_mat)
    dimnames(C) = list(dimnames(A)[[3]], dimnames(B)[[3]])
    C
}
arrayProd2tr <- function(A, B) {
    # p <- dim(A)[3]
    # k <- dim(B)[3]
    # # Preallocate output: k x n x n x p
    # result <- matrix(0,k, p,
    #                 dimnames = list(dimnames(B)[[3]], dimnames(A)[[3]]))
    # # Loop over p and k slices to avoid big intermediate arrays
    # for (ip in seq_len(p)) {
    #     for (ik in seq_len(k)) {
    #         result[ik, ip] <- tr(A[,,ip] %*% B[,,ik])
    #     }
    # }
    # result
    #High memory, but fast version
    p <- dim(A)[3]
    k <- dim(B)[3]
    n <- dim(A)[1]

    # reshape: each slice becomes a column
    Amat <- matrix(aperm(A, c(2,1,3)), n*n, p)  # vec(A[,,i]^T)
    Bmat <- matrix(B, n*n, k)                   # vec(B[,,j])

    # big matrix multiplication gives all traces at once
    result <- t(crossprod(Amat, Bmat)) # k x p

    dimnames(result) <- list(dimnames(B)[[3]], dimnames(A)[[3]])
    result
}
#Shift coordinates
shiftCoord = function(x, scaleByMax = TRUE){
    Mins = apply(vapply(x, FUN.VALUE = double(2), FUN = function(y){
        apply(y, 2, min)
    }), 1, min)
    Max = if(scaleByMax){
        max(vapply(x, FUN.VALUE = double(2), FUN = function(y){
             apply(y, 2, function(z) diff(range(z)))
        }))
    } else {1}
    lapply(x, function(y){
        y[, "x"] = (y[, "x"] - Mins[1])/Max
        y[, "y"] = (y[, "y"] - Mins[2])/Max
        y
    })
}
ptapply <- function(X, INDEX, FUN = NULL, ..., nCores = getOption("mc.cores", 2L)) {
    # Handle multiple grouping variables
    if (is.data.frame(INDEX) || is.list(INDEX)) {
        group_factor <- interaction(INDEX, drop = TRUE)
    } else {
        group_factor <- INDEX
    }
    # Split by group
    split_data <- split(X, group_factor)
    ngroups <- length(split_data)

    # Chunk groups across cores
    chunks <- split(seq_len(ngroups), cut(seq_len(ngroups), nCores, labels = FALSE))

    # Process chunks in parallel
    chunk_results <- parallel::mclapply(chunks, function(idx) {
        sapply(split_data[idx], FUN, ..., simplify = FALSE)
    }, mc.cores = nCores)
    out = unlist(chunk_results, recursive = FALSE)
    names(out) = names(split_data)
    return(out)
}
#' A wrapper for Matrix::bdiag maintaining names
#'
#' @param A,B Matrix to be used in \link[Matrix]{bdiag}
#' @return Same as \link[Matrix]{bdiag} but with dimnames
bdiagn = function(A, B){
    M <- bdiag(A, B)
    # Build new dimnames from components
    dimnames(M) <- list(c(rownames(A), rownames(B)), c(colnames(A), colnames(B)))
    M
}
getPower = function(x, sigLevel = 0.05){
    mean(x<sigLevel, na.rm  =TRUE)
}
#Fdrtol correction on p-values, return qvals
fdrtoolQval = function(pvals, verbose = FALSE, plot = FALSE,...){
    if(is.null(pvals)){
        warning("Null provided!")
        return(NULL)
    }
    tmp = fdrtool(c(pvals), statistic = "pvalue", plot = plot, verbose = verbose, ...)$qval
    if(is.matrix(pvals)){
        names(tmp) = sbivar:::makeNames(rownames(pvals), colnames(pvals))
    }
    tmp
}
#Extract correlated features (same names) for GP simulations
getSameNamesId = function(charVec, which = TRUE){
    tmp = gsub("x", "", gsub("y", "", gsub("X", "", gsub("Y", "", simplify2array(strsplit(charVec, "__"))))))
    out = tmp[1,]==tmp[2,]
    if(which)
        which(out)
    else
        out
}
getSensFDP = function(p, true){
    p[is.na(p)] = 1
    signif = p < sigLevel
    trueSignif = sum(signif & true, na.rm = TRUE)
    falseSignif = sum(signif & !true, na.rm = TRUE)
    FDP = if(any(signif>0)) sum(falseSignif)/sum(signif, na.rm = TRUE) else 0
    Sens = if(any(true)) sum(trueSignif)/sum(true) else NA
    c("FDP" = FDP, "Sens" = Sens)
}
getGenes = function(x){
   getSplitMat(x)[1,]
}
getMets = function(x){
    getSplitMat(x)[2,]
}
getSplitMat = function(x){
    sapply(x, function(y) strsplit(y, "__")[[1]])
}
scrambleRows = function(x){#Scramble matrix by rows, but retain rownames
    tmp = x[sample(nrow(x)),, drop = FALSE]
    rownames(tmp) = rownames(x)
    return(tmp)
}
normMat = function(x, norm, pseudoCount = 1e-8) {
    if (norm == "none") {
        out <- x
    } else {
        if (min(x) < 0) {
            warning("Normalization '", norm, "'is not recommended for real valued data!")
        }
        x <- x[id <- ((rs <- rowSums(x)) > 0), , drop = FALSE]
        dn <- dimnames(x)
        out <- if (norm == "rel") {
            x / rs[id]
        } else if (norm == "log") {
            log((x + pseudoCount) / rs[id])
        }
        dimnames(out) <- dn
    }
    colnames(out) <- make.names(colnames(out))
    return(out)
}
makeRowNames = function(x){
    apply(x[, c("Modality_X", "Modality_Y")], 1, paste, collapse = "__")
}
allNotSignif = function(x) {all(x>sigLevel, na.rm = TRUE)}
pasteCA = function(x, capital =FALSE){
    out = if(length(x)==1){
        paste(if(capital) "Sample" else "sample", x)
    } else {
        paste(if(capital) "Samples" else "samples", paste0(x[-length(x)], collapse = ", "), "and", tail(x, 1))
    }
    return(gsub("_", "\\_", out, fixed = TRUE))
}