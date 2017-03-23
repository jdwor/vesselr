#' @title 3D Volume Blobness
#' @description This function returns a blobness map for a 3D array or NIfTI volume. This blobness measure is based on the volume ratio described by Pierpaoli and Basser (1996).
#' @param image a 3D array or image of class \code{\link{nifti}}
#' @param color a string specifying whether blobs will appear darker ("dark") or brighter ("bright") than their surroundings
#' @param mask an array or \code{\link{nifti}} mask of voxels for which blobness will be calculated,
#' if \code{NULL} the blobness filter will be run for the full array.
#' Note that mask should be in the same space as the image volume
#' @param radius an integer specifying radius of the neighborhood (in voxels) for which the blobness should be calculated.
#' Note that this value essentially serves as the scale of the blob objects
#' @param parallel is a logical value that indicates whether the user's computer
#' is Linux or Unix (i.e. macOS), and should run the code in parallel
#' @return A 3D volume of the volume ratio blobness scores.
#' @examples \dontrun{
#' library(neurobase)
#' flair <- readnii('path/to/epi')
#' mask <- flair!=0
#' brightspots <- blobness3D(image = flair, mask = mask, radius = 5,
#'                       color = "bright", parallel = TRUE) }
#' @export
#' @references C. Pierpaoli, P.J. Basser (1996). Toward a Quantitative Assessment of Diffusion Anisotropy. Magnetic Resonance in Medicine. 36, pp. 893-906.
blobness3D=function(image, mask = NULL, radius = 5, color = "dark", parallel = FALSE){
  if(is.null(mask)){
    mask=image
    mask[mask==image]=1
  }

  eigvals=hessian3D(image,mask,radius,parallel)

  print("Calculating blobness measure")
  l1=eigvals$eigval1
  l2=eigvals$eigval2
  l3=eigvals$eigval3
  l1=as.vector(l1[mask==1])
  l2=as.vector(l2[mask==1])
  l3=as.vector(l3[mask==1])
  rm(eigvals)

  al1=abs(l1)
  al2=abs(l2)
  al3=abs(l3)

  bness=al1*al2*al3*((3/(al1+al2+al3))^3)

  if(color=="dark"){
    bness[l1<0 | l2<0 | l3<0] = 0
    bness[!is.finite(bness)] = 0
  }else if(color=="bright"){
    bness[l1>0 | l2>0 | l3>0] = 0
    bness[!is.finite(bness)] = 0
  }

  image[mask==1]<-bness
  return(image)
}
