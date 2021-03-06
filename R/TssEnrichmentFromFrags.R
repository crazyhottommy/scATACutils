#' checkClass function
#'
#' check whether the x object corresponds to the given class
#'
#' @param x object
#' @param class.name class name
#' @param var.name uses x object
#' @keywords internal
checkClass = function(x, class.name, var.name = deparse(substitute(x))){

  fun.name = match.call(call=sys.call(sys.parent(n=1)))[[1]]
  if(!class(x) %in% class.name)
    stop(paste(fun.name,': ',
               var.name,
               ' is not of class: ',
               paste(class.name, collapse=' '),
               '\n', sep=''))
}

### remove the tss that do not have coverage
### I took some code from the ScoreMatrix.R function in the genomation package.
### give the credit due :)
### see https://github.com/BIMSBbioinfo/genomation/blob/master/R/scoreMatrix.R#L113
constrainRanges = function(target, windows){

  checkClass(target, c('SimpleRleList','RleList','CompressedRleList'))
  checkClass(windows, 'GRanges')

  mcols(windows)$X_rank = 1:length(windows)
  r.chr.len = elementNROWS(target)
  constraint = GRanges(seqnames=names(r.chr.len),
                       IRanges(start=rep(1,length(r.chr.len)),
                               end=as.numeric(r.chr.len)))
  # suppressWarnings is done becuause GenomicRanges function give warnings
  #if you don't have the same seqnames in both objects
  win.list.chr = suppressWarnings(subsetByOverlaps(windows,
                                                   constraint,
                                                   type = "within",
                                                   ignore.strand = TRUE))

  if(length(win.list.chr) == 0)
    stop('All windows fell have coordinates outside windows boundaries')
  return(win.list.chr)
}



#' Calculate tss enrichment score from 10xscATAC fragments.tsv.gz file
#'
#' The memory usuage is not high (~5G), but takes ~2hours for 5000 cells.
#' @param frag_gz_file  fragment.tsv.gz file from 10x cellranger-atac output or
#' anyother tool but in the same format.
#' @param txs  a txdb object
#' @param flank flanking bp of tss (upstream and downstream)
#' @param endFlank  bp end flanks of flank for local noise control
#'
#'  -------flank---------flank-------
#'
#'  ----------------|-----------------
#'
#'  ---------------tss----------------
#'
#'  endFlank-----------------endFlank
#'
#'
#'
#' @param highest_tss_flank bp flanking tss windown for choosing the highest tss score.
#' The highest tss enrichment score is not always exactly at tss.
#' @param barcodeList A vector of valid barcodes
#' @param smooth window size to smooth
#' @param strand.aware consider tss strandness when calculating
#'
#' @return A dataframe with two columns: cells and the tss_score.
#' @export
#'
#' @examples
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(dplyr); library(readr); library(BiocParallel)
#' txs<- TxDb.Hsapiens.UCSC.hg19.knownGene
#' scores<- TssEnrichmentFromFrags("fragment.tsv.gz", txs = txs)
#' }

TssEnrichmentFromFrags <- function(frag_gz_file,
                                   txs,
                                   flank = 1000,
                                   endFlank = 100,
                                   highest_tss_flank= 50,
                                   smooth = 50,
                                   strand.aware = TRUE,
                                   workers = 1,
                                   barcodeList = NULL){

  # Make GRanges of fragments that are solid for the cells that we care about
  frags_valid <- data.table::fread(cmd = paste0("zcat < ", frag_gz_file)) %>%
    data.frame() %>%
    dplyr::mutate(V2 = V2 + 1) %>% # make it 1 based for R
    GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE)
  if (!is.null(barcodeList)){
    frags_valid<- frags_valid[frags_valid$V4 %in% barcodeList]
  }

  # common chromosome names, do it per cell instead, see TssEnrichmentSingleCell
  seqlev<- intersect(seqlevels(frags_valid), seqlevels(txs))
  frags_valid<- keepSeqlevels(frags_valid, seqlev, pruning.mode="coarse")

  # calculate coverage per cell
  frags_valid_per_cell<- split(frags_valid, frags_valid$V4)


  # this step can take minutes
  multicoreParam <- BiocParallel::MulticoreParam(workers = workers)
  # can add the chromosome length as additional argument for `coverage`
  # to get 0 coverages if there are no reads there.
  cvgs<- BiocParallel::bplapply(frags_valid_per_cell, function(x) coverage(x), BPPARAM = multicoreParam)
  txs<- transcripts(txs)
  txs <- unique(txs)

  txs.flanks<- promoters(txs, upstream = flank,
                         downstream = flank)
  txs.length<- length(txs.flanks)

  TssEnrichmentScores<- BiocParallel::bplapply(cvgs, TssEnrichmentSingleCell, txs.flanks, strand.aware = strand.aware, endFlank = endFlank, flank = flank, highest_tss_flank, smooth = smooth, BPPARAM = multicoreParam)

  enrichment<- do.call("rbind", TssEnrichmentScores)
  enrichment<- data.frame(cells = rownames(enrichment),
                          tss_score = enrichment[,1])
  return(enrichment)
}

TssEnrichmentSingleCell<- function(cvg, txs.flanks, strand.aware = TRUE, flank = 1000,
                                   endFlank = 100,
                                   highest_tss_flank= 50,
                                   smooth = 50 ){
  ## remove tss not in the coverage and assign a unique id for each tss: X_rank
  txs.flanks<- constrainRanges(cvg, txs.flanks)
  txs.length<- length(txs.flanks)
  if(length(txs.flanks)!=txs.length){
    warning(paste0(txs.length-length(txs.flanks),
                   " Tss removed because they fall out of the coverage"))
  }
  # common chromosomes
  chrs<- sort(intersect(names(cvg), as.character(unique(seqnames(txs.flanks)))))

  # convert GRanges to IntergerRangesList does not maintain the order
  # a unique id was given for each Ranges
  myViews<- Views(cvg[chrs],as(txs.flanks,"IntegerRangesList")[chrs]) # get subsets of RleList
  mat = lapply(myViews,function(x) t(viewApply(x,as.vector)) )
  mat = do.call("rbind",mat)

  r.list=split(mcols(txs.flanks)[,"X_rank"], as.vector(seqnames(txs.flanks))  )
  r.list=r.list[order(names(r.list))]
  ranks=do.call("c",r.list)
  rownames(mat) = ranks

  if(strand.aware == TRUE){
    orig.rows=txs.flanks[strand(txs.flanks) == '-',]$X_rank
    mat[rownames(mat) %in% orig.rows,] = mat[rownames(mat) %in%
                                               orig.rows, ncol(mat):1]
  }

  # reorder according to the original Granges (txs)
  mat = mat[order(ranks),]


  ### normlization by the endFlank local noise
  profile <- colSums(mat)
  profile_norm <- profile/mean(profile[c(1:endFlank,(flank*2-endFlank+1):(flank*2))])

  #smooth
  profile_norm_smooth <- zoo::rollmean(profile_norm, smooth, fill = 1)

  #enrichment
  max_finite <- function(x){
    suppressWarnings(max(x[is.finite(x)], na.rm=TRUE))
  }

  e <- max_finite(profile_norm_smooth[(flank-highest_tss_flank):(flank+highest_tss_flank)])
  return(e)
}
