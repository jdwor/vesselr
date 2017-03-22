#' @title NIfTI Gradient
#' @description This function returns the gradient images for a 3-dimensional NIfTI image.
#' @param image an image of class \code{\link{nifti}}
#' @param which a string specifying the gradient direction that should be returned. Either "all" for a list of x, y, and z gradients, or "x", "y", or "z" for a single image with the given gradient.
#'
#' @return Either a list of three gradient images or a single gradient image.
#' @examples \dontrun{
#' library(neurobase)
#' epi <- readnii('path/to/epi')
#' gradients <- niftiGradient(image = epi, which = "all") }
#' @export
#' @importFrom oro.nifti is.nifti
niftiGradient=function(image,which="all"){
  if(is.nifti(image)){
    if(which=="all"){
      dx=image
      dx@.Data[1,,]=0
      dx@.Data[dim(image)[1],,]=0
      dx@.Data[(2:(dim(image)[1]-1)),,]=(image@.Data[(3:dim(image)[1]),,]-
                                           image@.Data[(1:(dim(image)[1]-2)),,])/2

      dy=image
      dy@.Data[,1,]=0
      dy@.Data[,dim(image)[2],]=0
      dy@.Data[,(2:(dim(image)[2]-1)),]=(image@.Data[,(3:dim(image)[2]),]-
                                           image@.Data[,(1:(dim(image)[2]-2)),])/2

      dz=image
      dz@.Data[,,1]=0
      dz@.Data[,,dim(image)[3]]=0
      dz@.Data[,,(2:(dim(image)[3]-1))]=(image@.Data[,,(3:dim(image)[3])]-
                                           image@.Data[,,(1:(dim(image)[3]-2))])/2

      return(list(Dx=dx,Dy=dy,Dz=dz))
    }else if(which=="x"){
      dx=image
      dx@.Data[1,,]=0
      dx@.Data[dim(image)[1],,]=0
      dx@.Data[(2:(dim(image)[1]-1)),,]=(image@.Data[(3:dim(image)[1]),,]-
                                           image@.Data[(1:(dim(image)[1]-2)),,])/2
      return(dx)
    }else if(which=="y"){
      dy=image
      dy@.Data[,1,]=0
      dy@.Data[,dim(image)[2],]=0
      dy@.Data[,(2:(dim(image)[2]-1)),]=(image@.Data[,(3:dim(image)[2]),]-
                                           image@.Data[,(1:(dim(image)[2]-2)),])/2

      return(dy)
    }else if(which=="z"){
      dz=image
      dz@.Data[,,1]=0
      dz@.Data[,,dim(image)[3]]=0
      dz@.Data[,,(2:(dim(image)[3]-1))]=(image@.Data[,,(3:dim(image)[3])]-
                                           image@.Data[,,(1:(dim(image)[3]-2))])/2

      return(dz)
    }
  }else{
    print("Image must be a NifTI")
  }
}
