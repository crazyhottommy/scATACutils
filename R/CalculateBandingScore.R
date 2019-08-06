
## this is a re-write of https://github.com/shendurelab/mouse-atac/blob/master/banding_scores/calculate_nucleosome_banding_scores.R
## using tidyverse

#' Calculate banding score for a single cell
#'
#' @param df A dataframe with 3 columns. cell column contains the cell barcode for a single cell, insert_size column contains the insert size, read_count
#' column contains the number of paired reads for that insert size.
#'
#' @return A banding score for that cell
#' @export
#'
#' @examples
#'
get_banding_score<- function(df){
  df<- df %>%
    dplyr::filter(insert_size >0) %>%
    #make complete insert size from 1 to 1000
    tidyr::complete(insert_size = 1:1000, fill = list(read_count  = 0)) %>%
    dplyr::arrange(insert_size)
  periodogram<- spec.pgram(df$read_count / max(df$read_count), pad=0.3, tap=0.5, span=2, plot=F, fast=T)
  periodogram$freq<- 1/periodogram$freq

  banding_score<- sum(periodogram$spec[periodogram$freq >= 100 & periodogram$freq <= 300])
  return(banding_score)
}

#' Calculate Banding scores for all cells
#'
#' The input of this function is calculated by a python script:
#' \href{https://github.com/crazyhottommy/scATACtools/blob/master/python/get_insert_size_distribution_per_cell.py}{get_insert_size_distribution_per_cell.py} possorted_bam.bam pbmc_5k_insert_size.txt --barcodes barcodes.tsv
#'
#' @param df A dataframe with 3 columns. cell column contains the cell barcode for all cells, insert_size column contains the insert size, read_count
#' column contains the number of paired reads for that insert size
#' @param barcodeList  a vector of valid barcodes. e.g. the first column of barcodes.tsv from 10x
#'
#' @return Banding scores for all cells.
#' @export
#'
#' @examples
CalculateBandingScore<- function(df, barcodeList = NULL){
  ##  use furrr to speed up?

  if (!is.null(barcodeList)){
    banding_scores<- df %>%
      dplyr::filter(cell %in% barcodeList)
    dplyr::group_by(cell) %>%
      dplyr::nest() %>%
      dplyr::mutate(banding_score = map_dbl(data, get_banding_score)) %>%
      dplyr::select(-data)

  } else{
    banding_scores<- frags %>%
      dplyr::group_by(cell) %>%
      dplyr::nest() %>%
      dplyr::mutate(banding_score = map_dbl(data, get_banding_score)) %>%
      dplyr::select(-data)
  }

}





