# Create a correlation matrix by comparing 2 data matrices, or a single matrix to itself
# Arguments: xmatrix=first data matrix, 
#            ymatrix=second data matrix (optional, default=NULL),
#            display=will display correlation matrix heat map
#                    if set to TRUE (optional, default=TRUE),
#            xlabel=x axis label (optional, default="Axis 1")
#            ylabel=y axis label (optional, default="Axis 2")
#            cortype=correlation type, valid values are: "pearson","kendall" or "spearman" (optional, default="pearson")
GetCorrMatrix <- function(xmatrix, ymatrix=NULL, display=TRUE, xlabel="Axis 1", ylabel="Axis 2",cortype="pearson") {

  #Build correlation matrix from user defined matrices
  if(is.null(ymatrix)){ 

    #Compare a single matrix to itself
    COR <- cor(xmatrix,,,cortype)

    #Display heat map if display=true, otherwise just
    #return correlation matrix
    if(display) 
    {
       #Display heat map w/tick marks labeled as they are in the matrices
       image(x=seq(dim(xmatrix)[2]), y=seq(dim(xmatrix)[2]), z=COR, xlab=xlabel, ylab=ylabel,axes=FALSE)
       axis(1,at=seq(1,dim(xmatrix)[2],1),labels=colnames(xmatrix))
       axis(2,at=seq(1,dim(xmatrix)[2],1),labels=colnames(xmatrix))
       text(expand.grid(x=seq(dim(xmatrix)[2]), y=seq(dim(xmatrix)[2])), labels=round(c(COR),2))
    } 

  } else{

    #Compare two different matrices 	
    COR <- cor(xmatrix,ymatrix,,cortype)

    #Display heat map if display=true, otherwise just
    #return correlation matrix
    if(display) 
    {
       #Display heat map w/tick marks labeled as they are in the matrices
       image(x=seq(dim(xmatrix)[2]), y=seq(dim(ymatrix)[2]), z=COR, xlab=xlabel, ylab=ylabel,axes=FALSE)
       axis(1,at=seq(1,dim(xmatrix)[2],1),labels=colnames(xmatrix))
       axis(2,at=seq(1,dim(ymatrix)[2],1),labels=colnames(ymatrix))
       text(expand.grid(x=seq(dim(xmatrix)[2]), y=seq(dim(ymatrix)[2])), labels=round(c(COR),2))
    }

  }

  
  #Return correlation matrix
  return(COR)
}