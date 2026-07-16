factorizeDf = function(df){
    if(!is.null(df$Method)){
        df$Method = factor(df$Method, levels = methLevels, labels = methLabels, ordered = TRUE)
    }
    if(!is.null(df$spatPats)){
        df$spatPats = factor(df$spatPats, levels = spatPats, labels = spatPatsLabels, ordered = TRUE)
    }
    return(df)
}