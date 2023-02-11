library(tidyverse)
library(scp)
## Code adapted from Laurent Gato

## Read the data table
x <- read.delim("data/02ng/Task2-SearchTask/AllQuantifiedPeptides.tsv") |>
    select(-18) |> ## column 18 is all NAs
    janitor::clean_names() |>
    rename_with(~ gsub("intensity_ex_auto_", "int_", .x, fixed = TRUE)) |>
    rename_with(~ gsub("detection_type_ex_auto_", "det_", .x, fixed = TRUE))

colnames(x)[6] ="run1"
colnames(x)[7] ="run2"
colnames(x)[8] ="run3"
colnames(x)[9] ="run4"
colnames(x)[10] ="run5"
colnames(x)[11] ="run6"

## topN aggregating function
topN <- function(x, n = 3L) {
    ## Keep n rows with highest sum if there are > n rows available
    if (nrow(x) > n) {
        o <- order(rowSums(x, na.rm = TRUE),
                   decreasing = TRUE, na.last = TRUE)
        x <- x[head(o, n), ]
    }
    ## Aggregate by sum on top N features
    colSums(x, na.rm = TRUE)
}


## Create a QFeatures object, replace 0s by NA and aggregate peptides
## into proteins using sums, medians and sums on top 3.
qf <- readQFeatures(x, fnames = 1,
                     ecol = grep("^run", names(x)),
                     name = "peptides") |>
    zeroIsNA(i = "peptides") |>
    aggregateFeatures(i = "peptides",
                      fcol = "protein_groups",
                      name = "sums",
                      fun = colSums,
                      na.rm = TRUE) |>
    aggregateFeatures(i = "peptides",
                      fcol = "protein_groups",
                      name = "medians",
                      fun = colMedians,
                      na.rm = TRUE) |>
    aggregateFeatures(i = "peptides",
                      fcol = "protein_groups",
                      name = "top3",
                      fun = topN)


## Visualise peptide and protein intensities for three proteins of
## interest
qf |>
    filterFeatures(~!is.na(gene_names)) |>
    filterFeatures(~gene_names %in% c("PGM1", 'RRS1', 'ACO2')) |>
    longFormat(rowvars = "gene_names") |>
    as_tibble() |>
    ggplot(aes(x = colname, y = value, colour = rowname)) +
    geom_point() +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_line(aes(group = rowname)) +
    facet_grid(factor(gene_names, levels = c("ACO2", "RRS1", "PGM1")) ~
                   factor(assay, levels = c("peptides", "medians", "sums", "top3"))) +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(legend.position="none") +
    labs(y = "intensity", x = "quantitation method")

#note that there is no filtering in the peptide level
## sessionInfo()

## R version 4.2.2 Patched (2022-11-02 r83238)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Manjaro Linux

## Matrix products: default
## BLAS:   /usr/lib/libblas.so.3.11.0
## LAPACK: /usr/lib/liblapack.so.3.11.0

## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods
## [8] base

## other attached packages:
##  [1] scp_1.8.0                   QFeatures_1.9.0
##  [3] MultiAssayExperiment_1.24.0 SummarizedExperiment_1.28.0
##  [5] Biobase_2.58.0              GenomicRanges_1.50.1
##  [7] GenomeInfoDb_1.34.4         IRanges_2.32.0
##  [9] S4Vectors_0.36.1            BiocGenerics_0.44.0
## [11] MatrixGenerics_1.10.0       matrixStats_0.63.0
## [13] forcats_0.5.2               stringr_1.5.0
## [15] dplyr_1.0.10                purrr_0.3.5
## [17] readr_2.1.3                 tidyr_1.2.1
## [19] tibble_3.1.8                ggplot2_3.4.0
## [21] tidyverse_1.3.2

## loaded via a namespace (and not attached):
##  [1] ProtGenerics_1.30.0         bitops_1.0-7
##  [3] fs_1.5.2                    lubridate_1.9.0
##  [5] bit64_4.0.5                 httr_1.4.4
##  [7] tools_4.2.2                 backports_1.4.1
##  [9] utf8_1.2.2                  R6_2.5.1
## [11] DBI_1.1.3                   lazyeval_0.2.2
## [13] colorspace_2.0-3            withr_2.5.0
## [15] tidyselect_1.2.0            bit_4.0.5
## [17] compiler_4.2.2              textshaping_0.3.6
## [19] cli_3.4.1                   rvest_1.0.3
## [21] xml2_1.3.3                  DelayedArray_0.24.0
## [23] labeling_0.4.2              scales_1.2.1
## [25] systemfonts_1.0.4           XVector_0.38.0
## [27] pkgconfig_2.0.3             dbplyr_2.2.1
## [29] rlang_1.0.6                 readxl_1.4.1
## [31] generics_0.1.3              farver_2.1.1
## [33] jsonlite_1.8.4              vroom_1.6.0
## [35] googlesheets4_1.0.1         RCurl_1.98-1.9
## [37] magrittr_2.0.3              GenomeInfoDbData_1.2.9
## [39] Matrix_1.5-3                munsell_0.5.0
## [41] fansi_1.0.3                 MsCoreUtils_1.11.1
## [43] lifecycle_1.0.3             stringi_1.7.8
## [45] snakecase_0.11.0            MASS_7.3-58.1
## [47] zlibbioc_1.44.0             grid_4.2.2
## [49] parallel_4.2.2              crayon_1.5.2
## [51] lattice_0.20-45             haven_2.5.1
## [53] hms_1.1.2                   pillar_1.8.1
## [55] igraph_1.3.5                reprex_2.0.2
## [57] glue_1.6.2                  modelr_0.1.10
## [59] vctrs_0.5.1                 tzdb_0.3.0
## [61] cellranger_1.1.0            gtable_0.3.1
## [63] clue_0.3-63                 assertthat_0.2.1
## [65] BiocBaseUtils_1.0.0         janitor_2.1.0
## [67] broom_1.0.1                 AnnotationFilter_1.22.0
## [69] ragg_1.2.4                  googledrive_2.0.0
## [71] gargle_1.2.1                SingleCellExperiment_1.20.0
## [73] cluster_2.1.4               timechange_0.1.1
## [75] ellipsis_0.3.2