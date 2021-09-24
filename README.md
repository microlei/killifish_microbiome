Working title: Aromatic hydrocarbon pollution and its effect on Atlantic
killifish evolution and gut microbiota
================
Lei Ma
Last compiled on 23 September, 2021

## Abstract

The Atlantic killifish, Fundulus heteroclitus, is an abundant estuarine
fish broadly distributed along the eastern coast of the U.S. Known for
their hardiness to industrial pollution, populations of killifish have
repeatedly evolved tolerance to levels of aromatic hydrocarbon exposure
that would otherwise be lethal. This tolerance, which is heritable
across at least two generations, is linked to reduced activation of the
aryl hydrocarbon receptor (AHR) signaling pathway. In other animals such
as zebrafish and mice, the AHR has been shown to influence the gut
microbiome, particularly when activated by the model toxic pollutant
3,3’,4,4’,5-pendachlorobiphenyl (PCB-126). While there exists extensive
literature on killifish ecophysiology and evolution, little is known
about the natural microbiome of these fish, and there have been few
studies examining the downstream effects of their rapid adaptations on
their microbiota. In order to understand host and environmental effects
on killifish gut microbiota, we sampled two populations of wild fish -
New Bedford Harbor (NBH) and Scorton Creek (SC) - as well as lab reared
F2 generation fish originating from each of these wild populations. NBH
fish are known to have evolved tolerance to aromatic hydrocarbons while
SC fish are considered sensitive to aromatic hydrocarbon exposure. For
each fish, we used PCR to screen for common AHR-related genotypes and
used 16S rRNA-based amplicon sequencing to assess the microbial
composition of the gut microbiome. In preliminary results, we found that
fish from NBH, a more polluted site, had consistently higher microbial
alpha and beta diversity than fish from SC, a less polluted site.
Additionally, NBH fish guts were enriched in bacteria from the order
Vibrionales, which contain common fish pathogens. Overall, our results
represent a major first step in understanding the connection between the
evolutionary adaptation of killifish and their effects on the fish’s gut
microbiota.

## Directory description

### Data

Contains the metadata and data files that have been generated during
analysis phase. Ignores large files (like .rds files) used in some
scripts and notebooks due to github storage limits.

### Logs

Logs of the Snakemake dada2 data processing of the raw sequence reads

### Notebooks

Rmarkdown notebooks and generated during the exploratory phase of data
analysis. You can read the notebooks in your browser by clicking on the
.md files.

### Output

Output of the Snakemake dada2 data processing, including the ASV table
and some intermediate products. Does not include the quality profiles of
the samples.

### Scripts

R scripts used in the Snakemake workflow, in data analysis, and in
generating figures.

## Dependencies

For the dada2 Snakemake pipeline

-   See my [MiSeq](https://github.com/microlei/apprill-miseq) processing
    workflow for dependencies

For the data analysis and manuscript

-   See session info below

## Session Info

    ## R version 4.1.1 (2021-08-10)
    ## Platform: x86_64-apple-darwin17.0 (64-bit)
    ## Running under: macOS Mojave 10.14.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] skimr_2.1.3     decontam_1.12.0 vegan_2.5-7     lattice_0.20-44
    ##  [5] permute_0.9-5   phyloseq_1.36.0 forcats_0.5.1   stringr_1.4.0  
    ##  [9] dplyr_1.0.7     purrr_0.3.4     readr_2.0.1     tidyr_1.1.3    
    ## [13] tibble_3.1.4    ggplot2_3.3.5   tidyverse_1.3.1 here_1.0.1     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-153           bitops_1.0-7           fs_1.5.0              
    ##  [4] lubridate_1.7.10       httr_1.4.2             rprojroot_2.0.2       
    ##  [7] GenomeInfoDb_1.28.4    repr_1.1.3             tools_4.1.1           
    ## [10] backports_1.2.1        utf8_1.2.2             R6_2.5.1              
    ## [13] mgcv_1.8-36            DBI_1.1.1              BiocGenerics_0.38.0   
    ## [16] colorspace_2.0-2       rhdf5filters_1.4.0     ade4_1.7-18           
    ## [19] withr_2.4.2            tidyselect_1.1.1       compiler_4.1.1        
    ## [22] cli_3.0.1              rvest_1.0.1            Biobase_2.52.0        
    ## [25] xml2_1.3.2             scales_1.1.1           digest_0.6.27         
    ## [28] rmarkdown_2.11         XVector_0.32.0         base64enc_0.1-3       
    ## [31] pkgconfig_2.0.3        htmltools_0.5.2        dbplyr_2.1.1          
    ## [34] fastmap_1.1.0          rlang_0.4.11           readxl_1.3.1          
    ## [37] rstudioapi_0.13        generics_0.1.0         jsonlite_1.7.2        
    ## [40] RCurl_1.98-1.5         magrittr_2.0.1         GenomeInfoDbData_1.2.6
    ## [43] biomformat_1.20.0      Matrix_1.3-4           Rcpp_1.0.7            
    ## [46] munsell_0.5.0          S4Vectors_0.30.0       Rhdf5lib_1.14.2       
    ## [49] fansi_0.5.0            ape_5.5                lifecycle_1.0.0       
    ## [52] stringi_1.7.4          yaml_2.2.1             MASS_7.3-54           
    ## [55] zlibbioc_1.38.0        rhdf5_2.36.0           plyr_1.8.6            
    ## [58] grid_4.1.1             parallel_4.1.1         crayon_1.4.1          
    ## [61] splines_4.1.1          Biostrings_2.60.2      haven_2.4.3           
    ## [64] multtest_2.48.0        hms_1.1.0              knitr_1.34            
    ## [67] pillar_1.6.2           igraph_1.2.6           reshape2_1.4.4        
    ## [70] codetools_0.2-18       stats4_4.1.1           reprex_2.0.1          
    ## [73] glue_1.4.2             evaluate_0.14          data.table_1.14.0     
    ## [76] modelr_0.1.8           vctrs_0.3.8            tzdb_0.1.2            
    ## [79] foreach_1.5.1          cellranger_1.1.0       gtable_0.3.0          
    ## [82] assertthat_0.2.1       xfun_0.26              broom_0.7.9           
    ## [85] survival_3.2-13        iterators_1.0.13       IRanges_2.26.0        
    ## [88] cluster_2.1.2          ellipsis_0.3.2
