
#' Plot scatter plot of log10 depth vs FRIP or TSS score for scATACseq data
#'
#'
#' @param df A dataframe with at least 3 columns with names: cells, depth, FRIP or
#' cells, depth, tss_score
#' @param y x-axis will always be the log10(depth), y could be FRIP or tss_score
#' @param barcodes a vector of valid barcodes
#' @param hline dotted hline cutoff
#' @param vline dotted vline cutoff
#' @param col color of the dotted vline and hline
#'
#' @return a ggplot2 scatter plot
#' @export
#'
#' @examples
#' \dontrun{
#'
#' PlotScatter(frip, y = "FRIP")
#' PlotScatter(tss_scores, y = "tss_score")
#' }
PlotScatter<- function(df, y = "FRIP", barcodes = NULL,
                           hline = 0.6, vline = 3, col = "red"){
  if(!(all(c("cells", "depth") %in% colnames(df)) &&
       ("FRIP" %in% colnames(df) | "tss_score" %in% colnames(df)))) {
    stop("dataframe must contain column names: cells, depth, FRIP or
    cells, depth, tss_score.")
  }

  if (!is.null(barcodes)){
    if(y == "FRIP"){
      df<- df %>% dplyr::filter(FRIP <=1, cells %in% barcodes)

    } else if (y == "tss_score") {
      df <- df %>% dplyr::filter(cells %in% barcodes)
    } else {
      stop("only set y to FRIP or tss_score")
    }
  } else {
    if(y == "FRIP"){
      df<- df %>% dplyr::filter(FRIP <=1)
    } else if ( y == "tss_score"){
      invisible()
    } else {
      stop("only set y to FRIP or tss_score")
    }
  }

  p<- ggplot2::ggplot(df, ggplot2::aes(x= log10(depth), y = .data[[y]] )) +
    ggpointdensity::geom_pointdensity() +
    viridis::scale_color_viridis() +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(legend.position = "left") +
    ggplot2::geom_hline(yintercept = hline, linetype = 2, col = col, size = 1) +
    ggplot2::geom_vline(xintercept = vline, linetype = 2, col = col, size = 1) +
    ggplot2::xlab("log10(depth)") +
    ggplot2::ylab(y)

  p<- ggExtra::ggMarginal(p, type="density", fill = viridis::viridis(1))

  return(p)
}


