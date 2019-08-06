#' Extend a GRanges object upstream and downstream
#'
#' @param x A GRanges object
#' @param upstream bp for extending upstream
#' @param downstream bp for extending downstream
#'
#' @return An extended GRanges object
#' @export
#'
#' @examples
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


#' Add Gene sybmol to the GRanges object after calling genes(txdb)
#'
#' @param txdb A TxDb object, e.g. TxDb.Hsapiens.UCSC.hg19.knownGene
#' @param eg.db Either org.Hs.eg.db or org.Mm.eg.db
#'
#' @return An GRanges object for all the genes with gene symbol in the metadata column
#' @export
#'
#' @examples
addGeneNameToTxdb<- function(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                             eg.db = org.Hs.eg.db){
  gene<- GenomicFeatures::genes(txdb)
  ## 1: 1 mapping
  ss<- AnnotationDbi::select(eg.db, keys = gene$gene_id,
                             keytype="ENTREZID", columns = "SYMBOL" )
  gene$symbol<- ss[, 2]
  return(gene)
}


#' Plot scATACseq coverage for each cluster from tabix indexed fragment.tsv.gz file
#'
#' @param chrom chromosome
#' @param start chromosome start
#' @param end chromosome end
#' @param gene_name gene symbol for the gene one want's to plot. e.g. VEGFA. Use either chrom, start, end or
#' gene_name
#' @param upstream If use gene_name, specify the upstream bp for extending the GRanges
#' @param downstream If use gene_name, specify the downstream bp for extending the GRanges
#' @param fragment path to the fragment.tsv.gz file from 10x cellranger-atac output. indexed by tabix.
#' @param grouping path to a tsv file with three columns: "cell", "cluster", "depth". "cell" is the
#' cell barcode, "cluster" is the cluster id for each cell, and "depth" is the number of reads in that cell.
#' @param genome hg19, hg38 or mm9, mm10
#' @param txdb A TxDb object, e.g. TxDb.Hsapiens.UCSC.hg19.knownGene
#' @param eg.db Either org.Hs.eg.db or org.Mm.eg.db
#' @param ymax ymax for each track, if not specified, max of all the tracks will be calculated. Every
#' track will use the same ymax, so it is comparable across clusters.
#' @param label_cex size of the cluster label
#' @param yaxis_cex size of the y-axis
#' @param track_col color of the track
#' @param tick.dist chromosome tick distance to mark, default 10kb
#' @param minor.tick.dist minor chromosome tick distance to mark, default 2000 bp
#' @param tick_label_cex size of the tick label
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' library(org.Hs.eg.db)
#' library(org.Mm.eg.db)
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#' chrom<-  "chr12"
#' start<-  69730394
#' end<- 69760971
#' plotCoverageByGroup(chrom = chrom, start = start, end = end, fragment = "data/atac_viz/10k_pbmc/atac_v1_pbmc_10k_fragments.tsv.gz",
#'                     grouping = "data/atac_viz/grouping.txt", track_col = "red")

#' plotCoverageByGroup(gene_name = "NKG7", fragment = "data/atac_viz/10k_pbmc/atac_v1_pbmc_10k_fragments.tsv.gz",
#'                     grouping = "data/atac_viz/grouping.txt", tick_label_cex = 1, tick.dist = 5000,
#'                     minor.tick.dist = 1000)}

PlotCoverageByGroup<- function(chrom = NULL, start = NULL, end =NULL, gene_name, upstream = 2000,
                               downstream = 2000, fragment, grouping,
                               genome ='hg19', txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                               eg.db = org.Hs.eg.db,
                               ymax = NULL, label_cex = 1,
                               yaxis_cex = 1, track_col = "cadetblue2",
                               tick.dist = 10000, minor.tick.dist = 2000,
                               tick_label_cex = 1){
  grouping<- readr::read_tsv(grouping)
  if(! all(c("cell", "cluster", "depth") %in% colnames(grouping))) {
    stop('Grouping dataframe must have cell, cluster, and depth columns.')
  }
  ## get number of reads per group for normalization.
  ## not furthur normalize by the cell number in each group.
  grouping<-  grouping %>%
    group_by(cluster) %>%
    dplyr::mutate(cells_in_group = n(), total_depth_in_group = sum(depth)) %>%
    # reads per million (RPM)
    dplyr::mutate(scaling_factor = 1e6/(total_depth_in_group)) %>%
    ungroup() %>%
    dplyr::select(cell, cluster, scaling_factor)


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


  ## read in the fragment.tsv.gz file
  ## with "chr", "start", "end", "cell", "duplicate" columns. output from cellranger-atac
  # this returns a list
  reads<- Rsamtools::scanTabix(fragment, param = gr)

  reads<- reads[[1]] %>%
    tibble::enframe() %>%
    dplyr::select(-name) %>%
    tidyr::separate(value, into = c("chr", "start", "end", "cell", "duplicate"), sep = "\t") %>%
    dplyr::mutate_at(.vars = c("start", "end"), as.numeric) %>%
    # make it 1 based for R, the fragment.tsv is 0 based
    dplyr::mutate(start = start + 1) %>%
    inner_join(grouping) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

  ## plotting
  pp <- karyoploteR::getDefaultPlotParams(plot.type=1)
  pp$leftmargin <- 0.15
  pp$topmargin <- 15
  pp$bottommargin <- 15
  pp$ideogramheight <- 5
  pp$data1inmargin <- 10
  kp <- karyoploteR::plotKaryotype(genome = genome, zoom = gr, plot.params = pp)
  kp<- karyoploteR::kpAddBaseNumbers(kp, tick.dist = tick.dist, minor.tick.dist = minor.tick.dist,
                        add.units = TRUE, cex= tick_label_cex, digits = 6)
  ## calculate the normalized coverage
  normalized_coverage<- function(x){
    #if (!is(x, "GRangesList"))
    #  stop("'x' must be a GRangesList object")
    # specify the width to the whole chromosome to incldue the 0s
    cvgs<- lapply(x, function(x) coverage(x, width = kp$chromosome.lengths) * x$scaling_factor[1])
    return(cvgs)
  }
  # GRangesList object by group/cluster
  reads_by_group<- split(reads, reads$cluster)
  print(reads_by_group)
  coverage_norm<- normalized_coverage(reads_by_group)

  ## calculate the max coverage if not specified
  if (is.null(ymax)) {
    yaxis_common<- ceiling(max(lapply(coverage_norm, max) %>% unlist()))
  } else {
    yaxis_common<- ymax
  }
  ## add gene information
  genes.data <- karyoploteR::makeGenesDataFromTxDb(txdb,
                                      karyoplot=kp,
                                      plot.transcripts = TRUE,
                                      plot.transcripts.structure = TRUE)
  genes.data <- karyoploteR::addGeneNames(genes.data)
  genes.data <- karyoploteR::mergeTranscripts(genes.data)

  kp<- karyoploteR::kpPlotGenes(kp, data=genes.data, r0=0, r1=0.05, gene.name.cex = 1)

  for(i in seq_len(length(coverage_norm))) {
    read <- coverage_norm[[i]]
    at <- karyoploteR::autotrack(i, length(coverage_norm), r0=0.1, r1=1, margin = 0.1)
    karyoploteR::kpPlotCoverage(kp, data=read,
                         r0=at$r0, r1=at$r1, col = track_col, ymax = yaxis_common)
    karyoploteR::kpAxis(kp, ymin=0, ymax=yaxis_common, numticks = 2, r0=at$r0, r1=at$r1, cex = yaxis_cex, labels = c("", yaxis_common))
    karyoploteR::kpAddLabels(kp, labels = names(coverage_norm)[i], r0=at$r0, r1=at$r1,
                cex=label_cex, label.margin = 0.035)
  }
}
