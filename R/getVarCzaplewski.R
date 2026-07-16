getVarCzaplewski = function(x, y, prepCza){
    n = length(x)
    mxy = mean(x*y)
    mx2y2 = mean(x^2*y^2)
    varIxy = (mxy^2*n*(2*(1-prepCza["S2"] + prepCza["S1"]) +
                           (2*prepCza["S3"]-2*prepCza["S5"])*(n-3)+prepCza["S3"]*(n-2)*(n-3)) -
                  mx2y2 *(6*(1-prepCza["S2"]+prepCza["S1"]) +
                              (4*prepCza["S1"] - 2*prepCza["S2"])*(n-3) + prepCza["S1"]*(n-2)*(n-3)) +
                  n* ((1-prepCza["S2"]+prepCza["S1"])+(2*prepCza["S4"]-prepCza["S6"])*(n-3) + prepCza["S4"]*(n-2)*(n-3))
    )/((n-1)*(n-2)*(n-3)) - (mxy/(n-1))^2
    return(c("EIxy" = -mxy/(n-1), "varIxy" = unname(varIxy)))
}