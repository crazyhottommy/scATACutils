
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scATACutils: an R/Bioconductor package for working with 10x scATACseq data

<!-- badges: start -->

<!-- badges: end -->

The goal of scATACutils is to provide functions to work with 10x
scATACseq data. It includes functions to calculate qualtiy control
metrics such as banding score, Frip score, and TSS enrichment score etc.
Some convinient functions are provided for visulizations as well. Plot
the scatter plot of the Frip and read depth. Plot scATACseq coverage
tracks by group etc.

## Installation

You can install the released version of scATACutils from github with:

``` r
devtools::install_github("crazyhottommy/scATACutils")
```

## Example

Demonstration of some useful functions

``` r
library(scATACutils)
library(dplyr)
#> Warning: package 'dplyr' was built under R version 3.5.2
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(readr)
library(ggplot2)


## this takes 5 mins
#frip<- CalculateFripScore("~/5k_pbmc_atac/fragments.tsv.gz",
#                          "~/5k_pbmc_atac/peaks.bed")

frip<- read_tsv("~/5k_pbmc_atac/frip.txt", col_names = T)
#> Parsed with column specification:
#> cols(
#>   cells = col_character(),
#>   depth = col_double(),
#>   FRIP = col_double()
#> )

# a tibble with 3 columns, cell-barcode, depth and a Frip score
head(frip)
#> # A tibble: 6 x 3
#>   cells              depth  FRIP
#>   <chr>              <dbl> <dbl>
#> 1 AAACGAAAGAAAGCAG-1   283 0.113
#> 2 AAACGAAAGAAATACC-1     9 0.222
#> 3 AAACGAAAGAAATCTG-1     2 1    
#> 4 AAACGAAAGAACAGGA-1     5 0.8  
#> 5 AAACGAAAGAACCATA-1     1 1    
#> 6 AAACGAAAGAACGACC-1     3 0.333

## read in 5k pbmc atac data valid barcdoe
barcodes<- read_tsv("~/5k_pbmc_atac/pbmc_5k_atac_barcodes.tsv", col_names = F)
#> Parsed with column specification:
#> cols(
#>   X1 = col_character()
#> )

# the insert size distribution from https://github.com/crazyhottommy/scATACtools/blob/master/python/get_insert_size_distribution_per_cell.py

insert<- read_tsv("~/5k_pbmc_atac/pbmc_5k_insert_size.txt", col_names = T)
#> Parsed with column specification:
#> cols(
#>   cell = col_character(),
#>   insert_size = col_double(),
#>   read_count = col_double()
#> )

head(insert)
#> # A tibble: 6 x 3
#>   cell               insert_size read_count
#>   <chr>                    <dbl>      <dbl>
#> 1 TCACTCGAGTTCAAGA-1         231         82
#> 2 TCACTCGAGTTCAAGA-1          62        179
#> 3 TCACTCGAGTTCAAGA-1         273         51
#> 4 TCACTCGAGTTCAAGA-1          76        274
#> 5 TCACTCGAGTTCAAGA-1          56        228
#> 6 TCACTCGAGTTCAAGA-1         178        102

## this takes ~5mins
banding<- CalculateBandingScore(insert, barcodeList = NULL)

## distribution of the banding score after log10 transformation
ggplot(banding, aes(sample = log10(banding_score))) + 
  stat_qq() + 
  stat_qq_line(color = "red") +
  theme_bw(base_size = 14)
```

<img src="man/figures/README-example-1.png" width="60%" height="60%" />

#### TSS enrichment score

This can take ~2hours using 10 CPUs for 5000 cells.

``` r
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

tss_scores<- TssEnrichmentFromFrags("~/5k_pbmc_atac/fragments.tsv.gz",
                                    txs = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                    workers = 10,
                                    barcodeList = barcodes$X1)
```

### depth vs Frip and TSS score

read in a pre-computed tss score for the 5k pbmc atac dataset.

``` r
tss_scores<- readRDS("~/5k_pbmc_atac/5k_pbmc_atac_tss_scores.rds")
head(tss_scores)
#>                                 cells tss_score
#> AAACGAAAGCGCAATG-1 AAACGAAAGCGCAATG-1  8.915011
#> AAACGAAAGGGTATCG-1 AAACGAAAGGGTATCG-1  7.974064
#> AAACGAAAGTAACATG-1 AAACGAAAGTAACATG-1  6.430791
#> AAACGAAAGTTACACC-1 AAACGAAAGTTACACC-1  6.316003
#> AAACGAACAGAGATGC-1 AAACGAACAGAGATGC-1  8.408465
#> AAACGAACATGCTATG-1 AAACGAACATGCTATG-1  8.907742
head(frip)
#> # A tibble: 6 x 3
#>   cells              depth  FRIP
#>   <chr>              <dbl> <dbl>
#> 1 AAACGAAAGAAAGCAG-1   283 0.113
#> 2 AAACGAAAGAAATACC-1     9 0.222
#> 3 AAACGAAAGAAATCTG-1     2 1    
#> 4 AAACGAAAGAACAGGA-1     5 0.8  
#> 5 AAACGAAAGAACCATA-1     1 1    
#> 6 AAACGAAAGAACGACC-1     3 0.333

frip_tss<- inner_join(frip, tss_scores)
#> Joining, by = "cells"
#> Warning: Column `cells` joining character vector and factor, coercing into
#> character vector

head(frip_tss)
#> # A tibble: 6 x 4
#>   cells              depth  FRIP tss_score
#>   <chr>              <dbl> <dbl>     <dbl>
#> 1 AAACGAAAGCGCAATG-1 15317 0.834      8.92
#> 2 AAACGAAAGGGTATCG-1 17316 0.836      7.97
#> 3 AAACGAAAGTAACATG-1 25087 0.830      6.43
#> 4 AAACGAAAGTTACACC-1 22561 0.836      6.32
#> 5 AAACGAACAGAGATGC-1  8119 0.836      8.41
#> 6 AAACGAACATGCTATG-1  5261 0.735      8.91

PlotScatter(frip, y = "FRIP", barcodes = barcodes$X1, hline = 0.6,
            vline = 3)
```

<img src="man/figures/README-unnamed-chunk-3-1.png" width="60%" height="60%" />

TSS
score

``` r
PlotScatter(frip_tss, y = "tss_score", vline = 3, hline = 6)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="60%" height="60%" />

### Plot ATACseq tracks for each cluster of cells

``` r
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#> Loading required package: GenomicFeatures
#> Loading required package: BiocGenerics
#> Loading required package: parallel
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:parallel':
#> 
#>     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
#>     clusterExport, clusterMap, parApply, parCapply, parLapply,
#>     parLapplyLB, parRapply, parSapply, parSapplyLB
#> The following objects are masked from 'package:dplyr':
#> 
#>     combine, intersect, setdiff, union
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, append, as.data.frame, basename, cbind,
#>     colMeans, colnames, colSums, dirname, do.call, duplicated,
#>     eval, evalq, Filter, Find, get, grep, grepl, intersect,
#>     is.unsorted, lapply, lengths, Map, mapply, match, mget, order,
#>     paste, pmax, pmax.int, pmin, pmin.int, Position, rank, rbind,
#>     Reduce, rowMeans, rownames, rowSums, sapply, setdiff, sort,
#>     table, tapply, union, unique, unsplit, which, which.max,
#>     which.min
#> Loading required package: S4Vectors
#> Loading required package: stats4
#> 
#> Attaching package: 'S4Vectors'
#> The following objects are masked from 'package:dplyr':
#> 
#>     first, rename
#> The following object is masked from 'package:base':
#> 
#>     expand.grid
#> Loading required package: IRanges
#> 
#> Attaching package: 'IRanges'
#> The following objects are masked from 'package:dplyr':
#> 
#>     collapse, desc, slice
#> Loading required package: GenomeInfoDb
#> Loading required package: GenomicRanges
#> Loading required package: AnnotationDbi
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'AnnotationDbi'
#> The following object is masked from 'package:dplyr':
#> 
#>     select
library(org.Hs.eg.db)
#> 
PlotCoverageByGroup(gene_name = "MS4A1", fragment = "~/5k_pbmc_atac/atac_viz/10k_pbmc/atac_v1_pbmc_10k_fragments.tsv.gz",
                     grouping = "~/5k_pbmc_atac/atac_viz/grouping.txt", tick_label_cex = 1, tick.dist = 5000,
                     track_col = "red", 
                     label_cex = 0.5,
                     minor.tick.dist = 1000)
#> Parsed with column specification:
#> cols(
#>   cell = col_character(),
#>   cluster = col_character(),
#>   depth = col_double()
#> )
#> 'select()' returned 1:1 mapping between keys and columns
#> Joining, by = "cell"
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="60%" height="60%" />

## Acknowlegements

  - Thanks [Caleb](https://github.com/caleblareau) for sharing the FRIP
    score code.
  - Thanks [Ansu Satpathy](https://twitter.com/Satpathology) and
    [Jeffrey Granja](https://github.com/jeffmgranja) for sharing the TSS
    enrichment score codes. More details can be found at my blog post
    <https://divingintogeneticsandgenomics.rbind.io/post/calculate-scatacseq-tss-enrichment-score/>
    I referenced them in the source code.
  - The plotting track function is inspired by a post by Andrew Hill
    <http://andrewjohnhill.com/blog/2019/04/12/streamlining-scatac-seq-visualization-and-analysis/>
    and re-implemented using the
    [karyoploteR](http://bioconductor.org/packages/release/bioc/html/karyoploteR.html).
    Check [Signac](https://satijalab.org/signac/) by Tim Sturt for
    similar functionalities.
