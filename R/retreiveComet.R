#Do not build comet data, retreive them live
buildCoordsComet = function(fol, imgName = "0/0"){
    rawZar = zarr_open(paste0("Data/Duchini2026/UnumLocalia/UnumLocalia/", fol, "/comet/",
                              substr(fol, 1, 6), "_comet.ome.zarr"), mode = "r")
    img <- rawZar$get_item(imgName) #0/1-5 are downsamples, ignore these unless computationally necessary
    shape = dim(img$get_item(list(1, ":", 1, ":", ":"))$data)
    shape[4:5] = shape[4:5]*2^as.integer(substr(imgName, 3, 3)) #Correct for downsampling
    as.matrix(expand.grid("x" = seq_len(shape[4])-1, "y" = seq_len(shape[5])-1))
}
extractProteinComet = function(fol, protein, imgName = "0/0"){
    rawZar = zarr_open(paste0("Data/Duchini2026/UnumLocalia/UnumLocalia/", fol, "/comet/",
                              substr(fol, 1, 6), "_comet.ome.zarr"), mode = "r")
    img = rawZar$get_item(imgName)
    shape = dim(img$get_item(list(1, ":", 1, ":", ":"))$data)
    i = match(protein, proteinNames)-1
    out = c(img[1, i, 1, , ]$as.array())
    id = which(out > comet_thresholds[[fol]][make.names(protein), "Min"])
    list("out" = out[id], "id" = id)
}