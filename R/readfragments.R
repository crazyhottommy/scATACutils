
#' Read in the fragment.tsv.gz file into memory as a GenomicsRanges object
#'
#' @param frag_gz_file fragment.tsv.gz file from 10x cellranger output
#' @param cutSite whether use the cut sites instead of the whole fragmeng. Default TRUE for footprint
#'
#' @return A GenomicRanges object
#' @export
#'
#' @examples
#' \dontrun{
#' insertions<- ReadFragments("~/5k_pbmc_atac/fragments.tsv.gz", cutSite = TRUE)
#' }


ReadFragments<- function(frag_gz_file, cutSite = TRUE){
  fragments<- data.table::fread(cmd = paste0("zcat < ", frag_gz_file)) %>%
    data.frame() %>%
    dplyr::mutate(V2 = V2 + 1) %>%   # make it 1 based for R
    GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE)
  if (cutSite){
    insertions <- c(
      GenomicRanges::GRanges(seqnames = GenomicRanges::seqnames(fragments),
                             ranges = IRanges::IRanges(start(fragments), start(fragments)),
                             V4 = mcols(fragments)[, "V4"]),
      GenomicRanges::GRanges(seqnames = GenomicRanges::seqnames(fragments),
                             ranges = IRanges::IRanges(end(fragments), end(fragments)),
                             V4 = mcols(fragments)[, "V4"])
    )

  } else {
    insertions<-fragments
  }
  remove(fragments)
  gc()
  return(insertions)
}
