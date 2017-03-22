#' @title NIfTI Hessian
#' @description This function returns the eigenvalues of the hessian matrices for a 3-dimensional NIfTI volume.
#' @param image an image of class \code{\link{nifti}}
#' @param mask a mask of class \code{\link{nifti}},
#' if \code{NULL} the hessian filter will be run for the full array.
#' Note that mask should be in the same space as the image volume
#' @param parallel is a logical value that indicates whether the user's computer
#' is Linux or Unix (i.e. macOS), and should run the code in parallel
#' @return A list of three eigenvalue volumes.
#' @examples \dontrun{
#' library(neurobase)
#' epi <- readnii('path/to/epi')
#' mask <- epi!=0
#' hesseigs <- niftiHessian(image = epi, mask = mask) }
#' @export
#' @importFrom pbmcapply pbmclapply
#' @importFrom pbapply pblapply
#' @importFrom parallel detectCores
niftiHessian=function(image, mask = NULL, parallel = FALSE){
  if(is.null(mask)){
    mask=image
    mask[image]=1
  }

  print("Getting derivatives")
  grads=niftiGradient(image,which="all")
  gx=grads$Dx
  gy=grads$Dy
  gz=grads$Dz
  gradsx=niftiGradient(gx,which="all")
  gxx=gradsx$Dx
  gxy=gradsx$Dy
  gxz=gradsx$Dz
  gradsy=niftiGradient(gy,which="all")
  gyx=gradsy$Dx
  gyy=gradsy$Dy
  gyz=gradsy$Dz
  gradsz=niftiGradient(gz,which="all")
  gzx=gradsz$Dx
  gzy=gradsz$Dy
  gzz=gradsz$Dz

  rm(grads,gradsx,gradsy,gradsz,gx,gy,gz)

  print("Creating hessian matrices")
  bigmat=cbind(as.vector(gxx[mask==1]),as.vector(gxy[mask==1]),as.vector(gxz[mask==1]),
               as.vector(gyx[mask==1]),as.vector(gyy[mask==1]),as.vector(gyz[mask==1]),
               as.vector(gzx[mask==1]),as.vector(gzy[mask==1]),as.vector(gzz[mask==1]))

  rm(gxx,gxy,gxz,gyx,gyy,gyz,gzx,gzy,gzz)

  biglist=split(bigmat,row(bigmat))
  biglist=lapply(biglist,matrix,nrow=3,byrow=T)

  rm(bigmat)

  getevals=function(matrix){
    thiseig=eigen(matrix)$values
    sort=order(abs(thiseig))
    return(thiseig[sort])
  }

  print("Calculating eigenvalues")
  if(parallel==TRUE){
    result=matrix(unlist(pbmclapply(biglist,getevals,
                                    mc.cores=detectCores())),ncol=3,byrow=T)
  }else if(parallel==FALSE){
    result=matrix(unlist(pblapply(biglist,getevals)),ncol=3,byrow=T)
  }
  e1=mask
  e1[mask==1]<-result[,1]
  e2=mask
  e2[mask==1]<-result[,2]
  e3=mask
  e3[mask==1]<-result[,3]

  return(list(eigval1=e1,eigval2=e2,eigval3=e3))
}
