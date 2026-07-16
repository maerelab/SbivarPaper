filterNumber = function(x, minNum = 5){
    tab = table(x)
    names(tab)[tab >= minNum]
}
filterFeats = function(x, numFeats, include = NULL){
    out = if(numFeats >= ncol(x)){
        x
    } else {
        x[, rank(-colSums(x)) <= numFeats]
    }
    if(!is.null(include) && !all(include %in% colnames(out))){
        out = cbind(out, x[, intersect(colnames(x), setdiff(include, colnames(out)))])
    }
    return(out[, colSums(out)>0])
}
#Find most abundant features over a replicated dataset
filterFeatsList = function(List, numFeats, allPresent = FALSE, selectFrom = NULL){
    allFeats = unique(unlist(lapply(List, colnames)))
    mat = matrix(0, length(List), length(allFeats), dimnames = list(names(List), allFeats))
    for(sam in seq_along(List)){
        relExpr = colSums(List[[sam]])/sum(List[[sam]])
        mat[sam,names(relExpr)] = relExpr
    }
    vec = colMeans(mat)
    if(allPresent){
        vec = vec[colMeans(mat==0)==0]
    }
    if(!is.null(selectFrom)){
        vec = vec[selectFrom]
    }
    names(vec[rank(-vec) <= numFeats])
}
selMatching = function(mat, cn){
    mat[, intersect(cn, colnames(mat))]
}