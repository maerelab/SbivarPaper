#Wrap the CC.test
CCtestWrapper = function(x, y, xSeq, ySeq, tp, radius = .5,
                         correction = "variance", ...){
    imx = im(mat = x, xcol = xSeq, yrow = ySeq)
    imy = im(mat = y, xcol = xSeq, yrow = ySeq)
    CC.test(imx, imy, correction = correction, test.points = tp, radius = radius,
            verbose = FALSE, ...)
}