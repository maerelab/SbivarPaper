#Get the permutation p-value, remembering Phipson2010
getPermPval = function(stat, permStats){
    countPerm = sum(stat < permStats)
    nPerms <- length(permStats)
    if(countPerm > nPerms/2){
        countPerm = nPerms-countPerm
    }
    out = min((countPerm+1)/(nPerms+1)*2, 1)
    return(out)
}#Same, for a matrix and an array
getPermPvalMat = function(stat, permStats){
    countPerm = rowSums(c(stat) < permStats, dim = 2, na.rm = TRUE)
    nPerms <- dim(permStats)[3]
    id = countPerm > nPerms/2
    countPerm[id] = nPerms-countPerm[id]
    out = (countPerm+1)/(nPerms+1)*2
    out[out>1] = 1
    return(out)
}