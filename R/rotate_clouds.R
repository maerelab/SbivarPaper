rotate_clouds <- function(cloudList) {
    # Ensure matrices
    cloud1 <- as.matrix(cloudList[[1]])
    cloud2 <- as.matrix(cloudList[[2]])
    R <- get_rot_mat(cloud1)
    # Rotate both clouds
    cloud1_rot <- cloud1 %*% R
    cloud2_rot <- cloud2 %*% R
    colnames(cloud1_rot) = colnames(cloud2_rot) = c("x", "y")
    return(list(
        Stx = cloud1_rot,
        Msi = cloud2_rot
    ))
}
rotate_cloud <- function(cloud1) {
    cloud1_rot <- cloud1 %*% get_rot_mat(cloud1)
    colnames(cloud1_rot) = c("x", "y")
    cloud1_rot
}
get_rot_mat = function(x){
    cov_mat <- cov(x)
    # Eigen decomposition
    eig <- eigen(cov_mat)
    # Rotation matrix: eigenvectors
    # First eigenvector = direction of maximum variance (x-axis)
    eig$vectors
}