#Rotate a matrix, angle in radials
rotMat = function(angle){
    cbind(c(cos(angle), sin(angle)), c(-sin(angle), cos(angle)))
}
#Rotation
rotCoords = function(x, angle){
    x  %*% rotMat(angle)
}
#Mirroring
mirrorCoords = function(x){
    1-x
}
