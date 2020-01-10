addGeneNameToTxdb<- function(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                             eg.db = org.Hs.eg.db){
  gene<- GenomicFeatures::genes(txdb)
  ## 1: 1 mapping
  ss<- AnnotationDbi::select(eg.db, keys = gene$gene_id,
                             keytype="ENTREZID", columns = "SYMBOL" )
  gene$symbol<- ss[, 2]
  return(gene)
}

extend <- function(x, upstream=0, downstream=0)
{
  if (any(strand(x) == "*"))
    warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  trim(x)
}


#' Read in the fragment.tsv.gz file using Rsamtools::scanTabix for a specific region to
#' speed up
#'
#' @param the fragment.tsv.gz file need to be tabix indexed.
#' @param gr A GenomicRanges
#' @param barcodeList A vector of valid barcode to include
#' @param cutSite whether or not using the cutting sites instead of the whole fragment
#'
#' @return A GenomicRanges within the specified region
#' @export
#'
#' @examples
ReadFragmentsTabix<- function(fragment, gr, barcodeList, cutSite = FALSE){
  fragments<- Rsamtools::scanTabix(fragment, param = gr)
  fragments<- fragments[[1]] %>%
    tibble::enframe() %>%
    dplyr::select(-name) %>%
    tidyr::separate(value, into = c("chr", "start", "end", "cell", "duplicate"), sep = "\t") %>%
    filter(cell %in% barcodeList) %>%
    dplyr::mutate_at(.vars = c("start", "end"), as.numeric) %>%
    # make it 1 based for R, the fragment.tsv is 0 based
    dplyr::mutate(start = start + 1) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  if (cutSite){
    insertions <- c(
      GenomicRanges::GRanges(seqnames = GenomicRanges::seqnames(fragments),
                             ranges = IRanges::IRanges(start(fragments), start(fragments)),
                             cell = mcols(fragments)[, "cell"]),
      GenomicRanges::GRanges(seqnames = GenomicRanges::seqnames(fragments),
                             ranges = IRanges::IRanges(end(fragments), end(fragments)),
                             cell = mcols(fragments)[, "cell"])
    )

  } else {
    insertions<-fragments
  }
  remove(fragments)
  gc()
  return(insertions)
}



#' Plot raw count signal from scATACseq data per cell in heatmap
#'
#' @param chrom chromosome
#' @param start start
#' @param end end
#' @param gene_name the name of the gene to plot
#' @param upstream upstream bp to include
#' @param downstream downstream bp to include
#' @param fragment the fragment.tsv.gz file from 10x cellranger output
#' @param barcodeList a vector of valid cell barcode to include
#' @param genome hg19, hg38 or mm9, mm10
#' @param txdb a transcriptdb object. e.g.TxDb.Hsapiens.UCSC.hg19.knownGene
#' @param eg.db an orgdb object e.g.org.Hs.eg.db
#' @param cutSite whether use the cut sites instead of the whole fragmeng. Default FALSE
#' @param col_fun color of the corresponding count for c(0,1,2). e.g. c("white", "blue", "red")
#'
#' @return A complexHeatmap
#' @export
#'
#' @examples
#' \dontrun{
#' library(readr)
#' library(org.Hs.eg.db)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' barcodes<- read_tsv("~/5k_pbmc_atac/pbmc_5k_atac_barcodes.tsv", col_names = FALSE)
#' PlotCoverageByCell(gene_name = "MS4A1",
#' upstream = 2000,
#' downstream = 8000,
#' fragment= "~/5k_pbmc_atac/fragments.tsv.gz",
#' barcodeList=barcodes$X1,
#' genome = "hg19",
#' txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
#' eg.db = org.Hs.eg.db, cutSite = FALSE)
#' }
PlotCoverageByCell<- function(chrom = NULL, start = NULL, end = NULL, gene_name,
                              upstream = 2000, downstream = 2000, fragment, barcodeList,
                              genome = "hg19", txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                              eg.db = org.Hs.eg.db, cutSite = FALSE,
                              col_fun = c("white", "blue", "red")){

  if (is.null(chrom) & is.null(start) & is.null(end) & !is.null(gene_name)){
    gene <- addGeneNameToTxdb(txdb = txdb, eg.db = eg.db)
    gr<- gene[which(gene$symbol == gene_name)]
    if (length(gr) == 0){
      stop("gene name is not found in the database")
    } else if (length(gr) > 1) {
      gr<- gr[1]
      warning("multiple GRanges found for the gene, using the first one")
      gr<- extend(gr, upstream = upstream, downstream = downstream)
    } else {
      gr<- extend(gr, upstream = upstream, downstream = downstream)
    }
  } else if (!is.null(chrom) & !is.null(start) & !is.null(end)){
    gr<- GRanges(seq = chrom, IRanges(start = start, end = end ))
  }

  reads<- ReadFragmentsTabix(fragment = fragment, gr = gr,
                             cutSite = cutSite, barcodeList = barcodeList )

  reads_by_cell<- split(reads, reads$cell)

  #retain the seqlevel with only chromosomes in the GenomicRanges
  gr<- keepSeqlevels(gr, value=seqnames(gr)[1])
  chr_len<- seqlengths(gr)
  # add width as the length of the chromosome to add 0 coverge
  mats<- lapply(reads_by_cell, function(x) genomation::ScoreMatrixBin(target = GenomicRanges::coverage(x, width= chr_len), windows = gr, bin.num = width(gr)))
  col_fun<-  circlize::colorRamp2(c(0, 1, 2), col_fun)
  do.call("rbind", mats) %>%
    ComplexHeatmap::Heatmap(cluster_rows = FALSE,
                            cluster_columns = FALSE,
                            col = col_fun,
                            show_row_names = FALSE,
                            name = "count",
                            heatmap_legend_param = list(at=c(0,1,2), color_bar = "discrete", border= "black"))
}
