#' @title 3D Volume Hessian
#' @description This function returns the eigenvalues of the hessian matrices for a 3D array or NIfTI volume.
#' @param image a 3D array or image of class \code{\link{nifti}}
#' @param mask an array or \code{\link{nifti}} mask of voxels for which vesselness will be calculated,
#' with more selective masking improving speed significantly.
#' Note that mask should be in the same space as the image volume
#' @param radius an integer specifying radius of the neighborhood (in voxels) for which the hessian should be calculated
#' @param parallel is a logical value that indicates whether the user's computer
#' is Linux or Unix (i.e. macOS), and should run the code in parallel
#' @param cores if parallel = TRUE, cores is an integer value that indicates how many cores
#' the function should be run on
#'
#' @return A list of three eigenvalue volumes.
#' @examples \dontrun{
#' library(neurobase)
#' epi <- readnii('path/to/epi')
#' mask <- epi!=0
#' hesseigs <- hessian3D(image = epi, mask = mask) }
#' @export
#' @importFrom pbmcapply pbmclapply
#' @importFrom pbapply pblapply
#' @importFrom parallel detectCores
hessian3D = function (image, mask, radius = 1, parallel = FALSE, cores = 2) {
  print("Getting derivatives")
  # can avoid calculating some derivatives because of symmetrical Hessian
  grads = gradient3D(image, which = "all", radius = radius)
  gradsx = gradient3D(grads$Dx, which = "all", radius = radius)
  gxx = as.vector(gradsx$Dx[mask == 1])
  gxy = as.vector(gradsx$Dy[mask == 1])
  gxz = as.vector(gradsx$Dz[mask == 1])
  gyy = as.vector(gradient3D(grads$Dy, which = "y", radius = radius)[mask == 1])
  gyz = as.vector(gradient3D(grads$Dy, which = "z", radius = radius)[mask == 1])
  gzz = as.vector(gradient3D(grads$Dz, which = "z", radius = radius)[mask == 1])

  print("Creating hessian matrices")
  mask_voxels = sum(mask)
  # pre-allocating memory for bigmat
  bigmat = matrix(rep(0, 9 * mask_voxels),
                   nrow = mask_voxels, ncol = 9)
  bigmat[, 1] = gxx
  bigmat[, 2] = gxy
  bigmat[, 3] = gxz
  bigmat[, 4] = gxy
  bigmat[, 5] = gyy
  bigmat[, 6] = gyz
  bigmat[, 7] = gxz
  bigmat[, 8] = gyz
  bigmat[, 9] = gzz
  rm(grads, gradsx, gxx, gxy, gxz, gyy, gyz, gzz)

  # performs split() + lapply() in one step
  biglist = apply(bigmat, 1, matrix, nrow = 3, simplify = F)
  rm(bigmat)

  getevals = function(matrix) {
    thiseig = eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
    thiseig[order(abs(thiseig))]
  }

  print("Calculating eigenvalues")
  if (parallel == TRUE) {
    result = matrix(unlist(pbmclapply(biglist, getevals, mc.cores = cores)),
                    nrow = 3)
  }
  else if (parallel == FALSE) {
    result = matrix(unlist(pblapply(biglist, getevals)),
                    nrow = 3)
  }
  e1 = mask
  e1[mask == 1] = result[1, ]
  e2 = mask
  e2[mask == 1] = result[2, ]
  e3 = mask
  e3[mask == 1] = result[3, ]
  return(list(eigval1 = e1, eigval2 = e2, eigval3 = e3))
}
