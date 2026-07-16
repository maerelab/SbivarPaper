#Download a zip file, unzip it and remove the zip file
#Credit to https://stackoverflow.com/questions/3053833/using-r-to-download-zipped-data-file-extract-and-import-data
downloadAndUnzip = function(url, destfile, ...){
    tmp <- tempfile()
    download.file(url, destfile = tmp, ...)
    unzip(tmp, exdir = destfile)
    unlink(tmp)
}