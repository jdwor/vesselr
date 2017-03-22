#' @title NIfTI Vesselness
#' @description This function returns a vesselness image for a 3-dimensional NIfTI volume. This vesseless measure is based on the method described by Frangi (1998).
#' @param image an image of class \code{\link{nifti}}
#' @param vein.color a string specifying whether veins will appear darker ("dark") or brighter ("bright") than their surroundings
#' @param mask a mask of class \code{\link{nifti}},
#' if \code{NULL} the vesselness filter will be run for the full array.
#' Note that mask should be in the same space as the image volume
#' @return A nifti volume of the Frangi vesselness measure.
#' @examples \dontrun{
#' library(neurobase)
#' epi <- readnii('path/to/epi')
#' mask <- epi!=0
#' hesseigs <- niftiHessian(image = epi, vein.color = "dark", mask = mask) }
#' @export
#' @references A.F. Frangi, W.J. Niessen, K.L. Vincken, M.A. Viergever (1998). Multiscale vessel enhancement filtering. In Medical Image Computing and Computer-Assisted Intervention - MICCAI'98, W.M. Wells, A. Colchester and S.L. Delp (Eds.), Lecture Notes in Computer Science, vol. 1496 - Springer Verlag, Berlin, Germany, pp. 130-137.
niftiVesselness=function(image, vein.color = "dark", mask = NULL){
  if(is.null(mask)){
    mask=image
    mask[mask==image]=1
  }

  eigvals=niftiHessian(image,mask)

  print("Calculating vesselness measure")
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

  Ra=al2/al3
  Ra[!is.finite(Ra)]<-0
  Rb=al1/sqrt(al2*al3)
  Rb[!is.finite(Rb)]<-0

  S=sqrt(al1^2 + al2^2 + al3^2)
  A=2*(.5^2)
  B=2*(.5^2)
  C=2*(.5*max(S))^2

  rm(al1,al2,al3)

  eA=1-exp(-(Ra^2)/A)
  eB=exp(-(Rb^2)/B)
  eC=1-exp(-(S^2)/C)

  rm(Ra,Rb,S,A,B,C)

  vness=eA*eB*eC

  rm(eA,eB,eC)

  if(vein.color=="dark"){
    vness[l2<0 | l3<0] = 0
    vness[!is.finite(vness)] = 0
  }else if(vein.color=="bright"){
    vness[l2>0 | l3>0] = 0
    vness[!is.finite(vness)] = 0
  }

  image[mask==1]<-vness
  return(image)
}
