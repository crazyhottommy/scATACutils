

## this function copied and modified from https://github.com/caleblareau/scATAC_10X_raw_to_kmers/blob/master/example_kmers.R

#' Calculate Fragments Fraction of reads in peaks \href{https://www.encodeproject.org/data-standards/terms/#enrichment}{(FRiP)} score
#'
#' The fragments.tsv.gz file from 10x cellranger_atac output contains all barcodes found from
#' the bam file. one can provide a valid barcode vector to filter those out if needed.
#'
#' @param frag_gz_file  a fragment.tsv.gz file output from 10x cellranger-atac.
#' @param peaks a 3 column bed file for called peaks.
#'
#' @return a tibble with 3 columns: sample, depth, Frip score.
#' @export
#'
#' @examples
#'
#'\dontrun{
#' CalculateFripScore("fragments.tsv.gz", "peaks.bed")
#'}
CalculateFripScore<- function(frag_gz_file,
                              peaks, barcodeList = NULL){

  peaks_gr<- rtracklayer::import(peaks, format = "BED")

  if (is.null(barcodeList)){
    frags_valid <- data.table::fread(cmd = paste0("zcat < ", frag_gz_file)) %>%
      data.frame() %>%
      dplyr::mutate(V2 = V2 + 1) %>% # make it 1 based for R
      GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE)

  } else {
    frags_valid <- data.table::fread(paste0("zcat < ", frag_gz_file)) %>%
      data.frame() %>%
      dplyr::filter(V4 %in% barcodeList) %>%
      dplyr::mutate(V2 = V2 + 1) %>% # make it 1 based for R
      GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE)
  }


  # Get a denominator, per cell
  denom <- table(GenomicRanges::mcols(frags_valid)$V4)
  barcodes_found <- names(denom)

  # Get the overlaps with peaks
  ovPEAK <- GenomicRanges::findOverlaps(peaks_gr, frags_valid)

  # Establish a numeric index for the barcodes for sparse matrix purposes
  id <- factor(as.character(GenomicRanges::mcols(frags_valid)$V4), levels = barcodes_found)

  # Make sparse matrix with counts with peaks by  unique barcode
  countdf <- data.frame(peaks = S4Vectors::queryHits(ovPEAK),
                        sample = as.numeric(id)[S4Vectors::subjectHits(ovPEAK)]) %>%
    dplyr::group_by(peaks,sample) %>% dplyr::summarise(count = n()) %>% data.matrix()

  m <- Matrix::sparseMatrix(i = c(countdf[,1], length(peaks_gr)),
                            j = c(countdf[,2], length(barcodes_found)),
                            x = c(countdf[,3],0))
  colnames(m) <- barcodes_found

  # Make a polished colData
  colData <- tibble::tibble(
    cells = barcodes_found,
    depth = as.numeric(denom),
    FRIP = Matrix::colSums(m)/as.numeric(denom)
  )
  return(colData)
}
