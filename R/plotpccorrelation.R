
# This function is adapted and modified from a blog post by Andrew Hill
# http://andrewjohnhill.com/blog/2019/05/06/dimensionality-reduction-for-scatac-data/

#' Plot Correlation of PCs with a metadata column of a Seurat object or a SingleCellExperiment object
#'
#' @param object A Seurat or SingleCellExperiment Object
#' @param reduction The reduction name. e.g. pca, lsi, PCA, LSI
#' @param column the name of the metadata column or ColData
#' @param method correlation method, see cor
#'
#' @return A ggplot2 scatter plot
#' @export
#'
#' @examples
#' \dontrun{
#' PlotPCcorrelation(pbmc_seurat, reduction = "lsi")
#' }

PlotPCcorrelation<- function(object, reduction ="lsi", column = "nCount_ATAC",
                                                   method = "spearman"){
  if (class(object) == "Seurat"){
    stopifnot(reduction %in% names(object@reductions))
    stopifnot(column %in% colnames(object@meta.data))
    coords<- Seurat::Embeddings(object, reduction=reduction)
    column_value<-  object@meta.data[, column]
  } else if (class(object) == "SingleCellExperiment") {
    stopifnot(reduction %in% names(object@reducedDims))
    stopifnot(column %in% colnames(SingleCellExperiment::colData(object)))
    coords<- SingleCellExperiment::reducedDim(object, type = reduction)
    column_value<- SingleCellExperiment::colData(object)[, column]
  } else {
    stop("please only provide a Seurat object or SingleCellExperiment object")
  }

  correlations<-  apply(coords, 2, function(x) {cor(x, column_value, method= method)})
  # absolute value of correlation, the direction of PC is random
  correlations_df<- data.frame(correlation=abs(correlations), PC=1:ncol(coords))

  plot_obj<-  ggplot2::ggplot(correlations_df, ggplot2::aes(PC, correlation)) +
    ggplot2::geom_point() +
    ggplot2::theme_classic() +
    ggplot2::geom_hline(yintercept = 0, linetype='dashed', color='red')

  return(plot_obj)

}
