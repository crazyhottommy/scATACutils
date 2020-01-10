#' Plot Transcription factor motif footprint in heatmap and lineplot
#'
#' @param PWM TF position weight matrix
#' @param peaks A GenomicRanges for the ATACseq peaks
#' @param cvg coverage of the cutsite (from the fragment.tsv.gz file)
#' @param extend basepair to extend for the motif. 100bp by default
#' @param genome hg19, hg38 or mm9, mm10
#' @param heatmap_col heatmap color for min and max values
#' @param lineplot_cols lineplot color for the line and dotted vertical line
#'
#' @return a list of heatmap (EnrichedHeatmap object) and a lineplot (ggplot2 object)
#' @export
#'
#' @examples
#'
PlotMotifFootPrint<- function(PWM, peaks, cvg, extend = 100, genome= "hg19",
                              heatmap_col = c("white", "red"),
                              lineplot_cols = c("blue", "red")){
  motif_pos<- motifmatchr::matchMotifs(PWM, peaks, genome= genome, out = "position", bg= "genome")
  motif_len<- width(motif_pos[[1]])[1]
  mat<- genomation::ScoreMatrix(cvg, windows = extend(motif_pos[[1]], upstream = extend, downstream = extend), rpm = TRUE, strand.aware = TRUE)

  ## note this gives the same results but much slower, I will use ScoreMatrix instead
  ## mat<- EnrichedHeatmap::normalizeToMatrix(as(cvg, "GRanges"),
  ##                                 target = head(motif_pos[[1]], 8000),
  ##                                 extend = 20,
  ##                                 value_column = "score",
  ##                                 w = 1,
  ##                                 mean_mode = "absolute",
  ##                                 target_ratio = 0.065)

  mat<- EnrichedHeatmap::as.normalizedMatrix(t(scale(t(mat))), k_upstream = extend,
                                             k_downstream = extend,
                                             k_target = motif_len,
                                             extend = extend,
                                             signal_name = "atac",
                                             target_name = "motif")
  ## heatmap
  hp<- EnrichedHeatmap::EnrichedHeatmap(mat, name = "normalized\ncount",
                                        col = heatmap_col,
                                        use_raster =TRUE,
                                        raster_quality = 2)

  df<- colMeans(mat) %>%
    tibble::enframe(name = "pos", value = "cut_freq") %>%
    dplyr::mutate(pos = factor(pos, levels = pos))
  ## lineplot
  p<- ggplot2::ggplot(df, aes(x = pos, y = cut_freq)) +
    ggplot2::geom_line(group = 1, color = lineplot_cols[1]) +
    ggplot2::xlab(NULL) +
    ggplot2::theme_classic(base_size = 14) +
    ggplot2::scale_x_discrete(breaks = c("u1", "t1", paste0("t", motif_len), paste0("d", extend)),
                              labels =c (paste0(-extend, "bp"), "start", "end", paste0(extend, "bp"))) +
    ggplot2::geom_vline(xintercept = which(levels(df$pos) %in% c("t1", paste0("t", motif_len))), linetype = 2, color = lineplot_cols[2])

  return(list(heatmap = hp, lineplot= p))

}
