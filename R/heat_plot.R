#' Immune System Data
#
#'
#'
#' A data set used by De Vito et al. (2019) as an illustrative example.
#'
#' @format A list of two matrices, each with 63 columns and 285 and 140 rows, respectively.
#' @examples
#' \dontrun{
#' The commands below show how the dataset was obtained from libraries on the Bioconductor repository.
#' source("http://bioconductor.org/biocLite.R")
#' biocLite(c("limma", "curatedOvarianData", "RTCGAToolbox"), suppressUpdates=TRUE)
#' library(curatedOvarianData)
#' library(RTCGAToolbox)
#' data(package="curatedOvarianData")
#' data(GSE20565_eset)
#' data(GSE9891_eset)
#' im_response <- c(
#'  "TP53",
#'  "TUBA3C","PRKACG","FGF6","FGF23","FGF22","FGF20","ASB4","TUBB1","LAT","ULBP1","NCR1",
#'  "SIGLEC5","CD160","KLRD1","NCR3","TRIM9","FGF18","ICOSLG",
#'  "MYH2","C9","MBL2",
#'  "GRIN2B","POLR2F","CSF2","IL5","CRP","C8A","SPTA1","GRIN2A","CCR6","FGA","LBP",
#'  "DUSP9","FCN2","PRKCG","ADCY8","IL5RA","GRIN1","C8B",
#'  "GH2","TNFSF18","GRIN2D","FGB","PRL","SPTBN5","CD70","FGG","RASGRF1","IFNG",
#'  "SPTBN4","TRIM10","ACTN2","LTA","TNFSF11","GRIN2C","CAMK2B","SPTB","IL1A","TNFRSF13B",
#'  "ITGA2B","CAMK2A","TRIM31","EREG")
#' GSE20_eset <- t(as.matrix(GSE20565_eset))
#' GSE98_eset <- t(as.matrix(GSE9891_eset))
#' leD <- length(im_response)
#' mat1 <- matrix(sapply(1:leD, function(i) which(colnames(GSE20_eset)==im_response[i])))
#' mat2 <- matrix(sapply(1:leD, function(i) which(colnames(GSE98_eset)==im_response[i])))
#' GSE20 <- matrix(c(GSE20_eset[,mat1]), 140, leD)
#' GSE98 <-matrix(c(GSE98_eset[,mat2]), 285, leD)
#' data_immune <- list(GSE98, GSE20)}
#' "data_immune"
#' Generate heatmap of matrices
#' @param Matrix matrix to plot
#' @param limit limit for the values of the plot
#' @export
#' @examples
#' sigma = matrix(rbinom(25,1,.30), 5, 5)
#' plot.heat(sigma)
plot.heat <- function(Matrix,Xlab="",Ylab="",limit=c(-2,2)){
  Matrix = as.matrix(Matrix)
  colnames(Matrix)<-NULL
  rownames(Matrix)<-NULL
  x = reshape2::melt(data.frame(Matrix,ind=c(nrow(Matrix):1)),id="ind")
  colnames(x)<-c("X1","X2","value")
  p_X_heat = ggplot(data = x, aes(x=X2, y=X1, fill=value)) +
    theme_bw() +
    #geom_tile(show.legend = F) +
    geom_tile() +
    xlab(Xlab) +
    ylab(Ylab) +
    theme(axis.title=element_text(size=14,face="bold"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    scale_fill_gradient2(limits=limit) + 
    theme(legend.position="bottom")
  return(p_X_heat)
}