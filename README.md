Inflammation and Parkinson’s disease using Mendelian randomization
================

This repository contains the R scripts for reproducing the results of
the study “**A Mendelian randomization study investigating the causal
role of inflammation on Parkinson’s disease**”.

# Setup

Please install and download the latest R version (see
[here](https://www.r-project.org/)) and the latest RStudio version (see
[here](https://www.rstudio.com/products/rstudio/download/)).
Additionally, please install the latest version of the following R
packages on CRAN:

``` r
install.packages(
  c(
    "tidyverse", "readxl", "janitor", "rsnps", "LDlinkR", 
    "MendelianRandomization", "remotes"
  ), 
    dependencies = TRUE
)
```

The latest version of the following R packages on Bioconductor
repository:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
```

The latest version of the following R packages from source:

``` r
remotes::install_github(c(
  "MRCIEU/TwoSampleMR",
  "MRCIEU/MRInstruments",
  "gqi/MRMix",
  "rondolab/MR-PRESSO"
))
```

``` r
install_github(
  "WSpiller/MRPracticals",
  build_opts = c("--no-resave-data", "--no-manual"),
  build_vignettes = TRUE
)
```

The project contains one folder for each Mendelian randomization (MR)
analysis that was performed in the study. The following instructions can
be used to reproduce the results of each MR analysis. The sections will
be labelled as the folders used in the project.

Before running the scripts for reproducing the analysis, please make
sure to do the following:

1.  Open the RStudio project `mr.inflpd.Rproj`.

2.  Download GWAS PD data from
    [here](https://drive.google.com/file/d/1FZ9UL99LAqyWnyNBxxlx6qOUlfAnublN/view).
    Store the file in the folder `pd\data-raw` and then run the script
    `pd\data-raw\01_data_import.R` to import GWAS PD data.

3.  Download GWAS AAO data from
    [here](https://drive.google.com/file/d/1n-6eOF6gaIxP9dLHx_QCndeaOQU2uhCf/view).
    Store the file in the folder `aao\data-raw`, unzip the file and run
    the script `aao\data-raw\01_data_import.R` to import GWAS AAO data.

4.  Set up a Personal Access Token to access LDlink API via `LDlinkR` R
    package. Instructions on how to do it can be found
    [here](https://cran.r-project.org/web/packages/LDlinkR/vignettes/LDlinkR.html).

# ahluwalia_il6

The folder contains the R scripts that can be used to reproduce the MR
analysis on IL-6 using data from [Ahluwalia et
al. (2021)](https://academic.oup.com/hmg/article/30/5/393/6124523)
study.

Please, run the scripts as follows:

1.  Download the supplementary tables file from
    [here](https://academic.oup.com/hmg/article/30/5/393/6124523#supplementary-data)
    (IL6_GWAS_Supp_Tables_HMG_11012021_ddab023 - xlsx file). Then copy
    the data in the sheet labelled as “Table 1. IL6 GWAS sig. SNPs” and
    store them into `ahluwalia_il6\data-raw\instruments_ahluwalia.xlsx`

2.  Run the script `ahluwalia_il6\data-raw\01_import_instruments.R` to
    import the genetic associations data.

3.  Run the script `ahluwalia_il6\data\02_tidy_instruments.R` to clean
    and select genetic variants for the analysis.

4.  Run the script `ahluwalia_il6\data\03_align_aoo_effects.R` to
    harmonize and prepare IL-6 and AAO data for MR analysis.

5.  Run the script `ahluwalia_il6\data\03_align_PD_effects.R` to
    harmonize and prepare IL-6 and PD data for MR analysis.

6.  Run the script `ahluwalia_il6\analysis\04_il6_aoo.R` to perform MR
    analyses for IL-6 and AAO.

7.  Run the script `ahluwalia_il6\analysis\04_il6_pd.R` to perform MR
    analyses for IL-6 and PD.

# il6r_il6

The folder contains the R scripts that can be used to reproduce the MR
analysis on IL-6 using data from [IL6R Genetics Consortium and Emerging
Risk Factors Collaboration
(2021)](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(11)61931-4/fulltext)
study.

Please, run the scripts as follows:

1.  Store the genetic association data of `rs2228145` with IL6 into
    `il6r_il6\data-raw\instruments_il6r.xlsx`

2.  Run the script `il6r_il6\data-raw\01_import_instruments.R` to import
    the genetic associations data.

3.  Run the script `il6r_il6\data\02_tidy_instruments.R` to clean and
    select genetic variants for the analysis.

4.  Run the script `il6r_il6\data\03_align_aoo_effects.R` to harmonize
    and prepare IL-6 and AAO data for MR analysis.

5.  Run the script `il6r_il6\data\03_align_PD_effects.R` to harmonize
    and prepare IL-6 and PD data for MR analysis.

6.  Run the script `il6r_il6\analysis\04_il6_aoo.R` to perform MR
    analyses for IL-6 and AAO.

7.  Run the script `il6r_il6\analysis\04_il6_pd.R` to perform MR
    analyses for IL-6 and PD.

# folkersen_il6

The folder contains the R scripts that can be used to reproduce the MR
analysis on IL-6 using data from [Folkersen et
al. (2020)](https://www.nature.com/articles/s42255-020-00287-2) study.

Please, run the scripts as follows:

1.  Summary-statistics from Folkersen et al. GWAS can be found here
    [here](https://zenodo.org/record/2615265#.Ycrj71nSKUk). Download the
    `IL-6.txt.gz` file and store it into the folder
    `folkersen_il6\data-raw`

2.  Run the script `folkersen_il6\data-raw\01_data_import.R` to import
    the genetic associations data.

3.  Run the script `folkersen_il6\data-raw\02_tidy_instruments.R` to
    tidy summary-level GWAS data.

4.  Run the script `folkersen_il6\data\03_prepare_instruments.R` to
    instruments data for MR analysis.

5.  Run the script `folkersen_il6\data\04_align_PD_effects.R` to
    harmonize and prepare IL-6 and PD data for MR analysis.

6.  Run the script `folkersen_il6\data\04_align_aao_effects.R` to
    harmonize and prepare IL-6 and AAO data for MR analysis.

7.  Run the script `folkersen_il6\analysis\05_il6_aoo.R` to perform MR
    analyses for IL-6 and AAO.

8.  Run the script `folkersen_il6\analysis\05_il6_pd.R` to perform MR
    analyses for IL-6 and PD.

# ccgc_crp

The folder contains the R scripts that can be used to reproduce the MR
analysis on CRP using data from [CCGC
(2011)](https://www.bmj.com/content/342/bmj.d548) study.

Please, run the scripts as follows:

1.  Store the genetic association data from Figure 2 in the manuscript
    into `ccgc_crp\data-raw\instruments_ccgc.xlsx`

2.  Run the script `ccgc_crp\data-raw\01_import_instruments.R` to import
    the genetic associations data.

3.  Run the script `ccgc_crp\data\02_tidy_instruments.R` to clean and
    select genetic variants for the analysis.

4.  Run the script `ccgc_crp\data\03_align_aoo_effects.R` to harmonize
    and prepare CRP and AAO data for MR analysis.

5.  Run the script `ccgc_crp\data\03_align_pd_effects.R` to harmonize
    and prepare CRP and PD data for MR analysis.

6.  Run the script `ccgc_crp\analysis\04_il6_aoo.R` to perform MR
    analyses for CRP and AAO.

7.  Run the script `ccgc_crp\analysis\04_il6_PD.R` to perform MR
    analyses for CRP and PD.

# ligthart_crp

The folder contains the R scripts that can be used to reproduce the MR
analysis on CRP using data from [Ligthart et
al. (2018)](https://www.sciencedirect.com/science/article/pii/S0002929718303203)
study.

Please, run the scripts as follows:

1.  Store the genetic association data from the manuscript (Table 1) and
    Supplementary Table S3 into
    `ligthart_crp\data-raw\instruments_Ligthart.xlsx`

2.  Run the script `ligthart_crp\data-raw\01_import_instruments.R` to
    import the genetic associations data.

3.  Run the script `ligthart_crp\data\02_tidy_instruments.R` to clean
    and select genetic variants for the analysis.

4.  Run the script `ligthart_crp\data\03_align_aoo_effects.R` to
    harmonize and prepare CRP and AAO data for MR analysis.

5.  Run the script `ligthart_crp\data\03_align_PD_effects.R` to
    harmonize and prepare CRP and PD data for MR analysis.

6.  Run the script `ligthart_crp\analysis\04_il6_aoo.R` to perform MR
    analyses for CRP and AAO.

7.  Run the script `ligthart_crp\analysis\04_il6_pd.R` to perform MR
    analyses for CRP and PD.

# folkersen_il1ra

The folder contains the R scripts that can be used to reproduce the MR
analysis on IL-6 using data from [Folkersen et
al. (2020)](https://www.nature.com/articles/s42255-020-00287-2) study.

Please, run the scripts as follows:

1.  Summary-statistics from Folkersen et al. GWAS can be found here
    [here](https://zenodo.org/record/2615265#.Ycrj71nSKUk). Download the
    `IL-1ra.txt.gz` file and store it into the folder
    `folkersen_il1ra\data-raw`

2.  Run the script `folkersen_il1ra\data-raw\01_data_import.R` to import
    the genetic associations data.

3.  Run the script `folkersen_il1ra\data-raw\02_tidy_instruments.R` to
    tidy summary-level GWAS data.

4.  Run the script `folkersen_il1ra\data\03_prepare_instruments.R` to
    instruments data for MR analysis.

5.  Run the script `folkersen_il1ra\data\04_align_PD_effects.R` to
    harmonize and prepare IL-1ra and PD data for MR analysis.

6.  Run the script `folkersen_il1ra\data\04_align_aao_effects.R` to
    harmonize and prepare IL-1ra and AAO data for MR analysis.

7.  Run the script `folkersen_il1ra\analysis\05_il1ra_aoo.R` to perform
    MR analyses for IL-1ra and AAO.

8.  Run the script `folkersen_il1ra\analysis\05_il1ra_pd.R` to perform
    MR analyses for IL-1ra and PD.

# herder_il1ra

The folder contains the R scripts that can be used to reproduce the MR
analysis on IL-1ra using data from [Herder et
al. (2014)](https://diabetesjournals.org/diabetes/article/63/12/4343/40454/Genetic-Determinants-of-Circulating-Interleukin-1)
study.

Please, run the scripts as follows:

1.  Store the genetic association data from the manuscript (Table 2)
    into `herder_il1ra\data-raw\instruments_herder.xlsx`

2.  Run the script `herder_il1ra\data-raw\01_import_instruments.R` to
    import the genetic associations data.

3.  Run the script `herder_il1ra\data\02_tidy_instruments.R` to clean
    and select genetic variants for the analysis.

4.  Run the script `herder_il1ra\data\03_align_aoo_effects.R` to
    harmonize and prepare IL-1ra and AAO data for MR analysis.

5.  Run the script `herder_il1ra\data\03_align_PD_effects.R` to
    harmonize and prepare IL-1ra and PD data for MR analysis.

6.  Run the script `herder_il1ra\analysis\04_il1ra_aoo.R` to perform MR
    analyses for IL-1ra and AAO.

7.  Run the script `il6r_il6\analysis\04_il1ra_pd.R` to perform MR
    analyses for IL-1ra and PD.

# prins_tnfalpha

The folder contains the R scripts that can be used to reproduce the MR
analysis on TNF-alpha using data from [Prins
(2016)](https://research.rug.nl/en/publications/inflammatory-biomarker-genomics-from-discovery-to-causality)
study.

Please, run the scripts as follows:

1.  Store the genetic association data from Table 2 in the manuscript
    into `prins_tnfalpha\data-raw\instruments_tnf.xlsx`. The manuscript
    can be found
    [here](https://pure.rug.nl/ws/portalfiles/portal/35098958/Chapter_3_.pdf).

2.  Run the script `prins_tnfalpha\data-raw\01_import_instruments.R` to
    import the genetic associations data.

3.  Run the script `prins_tnfalpha\data\02_tidy_instruments.R` to clean
    and select genetic variants for the analysis.

4.  Run the script `prins_tnfalpha\data\03_align_aoo_effects.R` to
    harmonize and prepare TNF-alpha and AAO data for MR analysis.

5.  Run the script `prins_tnfalpha\data\03_align_pd_effects.R` to
    harmonize and prepare TNF-alpha and PD data for MR analysis.

6.  Run the script `prins_tnfalpha\analysis\04_tnf_aoo.R` to perform MR
    analyses for TNF-alpha and AAO.

7.  Run the script `prins_tnfalpha\analysis\04_tnf_PD.R` to perform MR
    analyses for TNF-alpha and PD.

# Workflow reproducibility

Reproducibility of the analyses can be achieved using R version 4.1.2
and the following packages versions:

``` r
library(tidyverse)
#> -- Attaching packages --------------------------------------- tidyverse 1.3.1 --
#> v ggplot2 3.3.5     v purrr   0.3.4
#> v tibble  3.1.6     v dplyr   1.0.8
#> v tidyr   1.2.0     v stringr 1.4.0
#> v readr   2.1.2     v forcats 0.5.1
#> -- Conflicts ------------------------------------------ tidyverse_conflicts() --
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
library(MendelianRandomization)
#> Warning in .recacheSubclasses(def@className, def, env): undefined subclass
#> "numericVector" of class "Mnumeric"; definition not updated
library(TwoSampleMR)
#> TwoSampleMR version 0.5.6 
#> [>] New: Option to use non-European LD reference panels for clumping etc
#> [>] Some studies temporarily quarantined to verify effect allele
#> [>] See news(package='TwoSampleMR') and https://gwas.mrcieu.ac.uk for further details
#> 
#> Attaching package: 'TwoSampleMR'
#> The following objects are masked from 'package:MendelianRandomization':
#> 
#>     mr_ivw, mr_median
library(knitr)
library(meta)
#> Loading 'meta' package (version 5.2-0).
#> Type 'help(meta)' for a brief overview.
#> Readers of 'Meta-Analysis with R (Use R!)' should install
#> older version of 'meta' package: https://tinyurl.com/dt4y5drs
library(officer)
library(flextable)
#> Warning: package 'flextable' was built under R version 4.1.3
#> 
#> Attaching package: 'flextable'
#> The following object is masked from 'package:purrr':
#> 
#>     compose
library(patchwork)
library(readxl)
#> 
#> Attaching package: 'readxl'
#> The following object is masked from 'package:officer':
#> 
#>     read_xlsx
library(janitor)
#> 
#> Attaching package: 'janitor'
#> The following objects are masked from 'package:stats':
#> 
#>     chisq.test, fisher.test
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
#> Loading required package: BSgenome
#> Loading required package: BiocGenerics
#> 
#> Attaching package: 'BiocGenerics'
#> The following object is masked from 'package:flextable':
#> 
#>     width
#> The following objects are masked from 'package:dplyr':
#> 
#>     combine, intersect, setdiff, union
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
#>     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
#>     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
#>     union, unique, unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> Loading required package: stats4
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:MendelianRandomization':
#> 
#>     values
#> The following objects are masked from 'package:dplyr':
#> 
#>     first, rename
#> The following object is masked from 'package:tidyr':
#> 
#>     expand
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> 
#> Attaching package: 'IRanges'
#> The following object is masked from 'package:TwoSampleMR':
#> 
#>     trim
#> The following objects are masked from 'package:dplyr':
#> 
#>     collapse, desc, slice
#> The following object is masked from 'package:purrr':
#> 
#>     reduce
#> The following object is masked from 'package:grDevices':
#> 
#>     windows
#> Loading required package: GenomeInfoDb
#> Loading required package: GenomicRanges
#> Loading required package: Biostrings
#> Loading required package: XVector
#> 
#> Attaching package: 'XVector'
#> The following object is masked from 'package:purrr':
#> 
#>     compact
#> 
#> Attaching package: 'Biostrings'
#> The following object is masked from 'package:base':
#> 
#>     strsplit
#> Loading required package: rtracklayer
library(rsnps)
library(openxlsx)
library(LDlinkR)
library(MRMix)
library(MRPRESSO)
library(MRPracticals)

si <- sessionInfo()
print(si, locale = FALSE)
#> R version 4.1.2 (2021-11-01)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 10 x64 (build 19044)
#> 
#> Matrix products: default
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] MRPracticals_0.0.1                      
#>  [2] MRPRESSO_1.0                            
#>  [3] MRMix_0.1.0                             
#>  [4] LDlinkR_1.1.2                           
#>  [5] openxlsx_4.2.5                          
#>  [6] rsnps_0.5.0.0                           
#>  [7] SNPlocs.Hsapiens.dbSNP144.GRCh37_0.99.20
#>  [8] BSgenome_1.62.0                         
#>  [9] rtracklayer_1.54.0                      
#> [10] Biostrings_2.62.0                       
#> [11] XVector_0.34.0                          
#> [12] GenomicRanges_1.46.0                    
#> [13] GenomeInfoDb_1.30.0                     
#> [14] IRanges_2.28.0                          
#> [15] S4Vectors_0.32.2                        
#> [16] BiocGenerics_0.40.0                     
#> [17] janitor_2.1.0                           
#> [18] readxl_1.3.1                            
#> [19] patchwork_1.1.1                         
#> [20] flextable_0.7.0                         
#> [21] officer_0.4.1                           
#> [22] meta_5.2-0                              
#> [23] knitr_1.37                              
#> [24] TwoSampleMR_0.5.6                       
#> [25] MendelianRandomization_0.6.0            
#> [26] forcats_0.5.1                           
#> [27] stringr_1.4.0                           
#> [28] dplyr_1.0.8                             
#> [29] purrr_0.3.4                             
#> [30] readr_2.1.2                             
#> [31] tidyr_1.2.0                             
#> [32] tibble_3.1.6                            
#> [33] ggplot2_3.3.5                           
#> [34] tidyverse_1.3.1                         
#> 
#> loaded via a namespace (and not attached):
#>   [1] uuid_1.0-4                  backports_1.4.1            
#>   [3] systemfonts_1.0.4           plyr_1.8.6                 
#>   [5] lazyeval_0.2.2              splines_4.1.2              
#>   [7] mr.raps_0.2                 gmp_0.6-4                  
#>   [9] BiocParallel_1.28.0         digest_0.6.29              
#>  [11] foreach_1.5.2               htmltools_0.5.2            
#>  [13] fansi_1.0.2                 magrittr_2.0.2             
#>  [15] tzdb_0.2.0                  modelr_0.1.8               
#>  [17] matrixStats_0.61.0          colorspace_2.0-3           
#>  [19] rvest_1.0.2                 haven_2.4.3                
#>  [21] xfun_0.30                   crayon_1.5.0               
#>  [23] RCurl_1.98-1.6              jsonlite_1.8.0             
#>  [25] lme4_1.1-28                 survival_3.2-13            
#>  [27] iterators_1.0.14            glue_1.6.2                 
#>  [29] gtable_0.3.0                zlibbioc_1.40.0            
#>  [31] MatrixModels_0.5-0          DelayedArray_0.20.0        
#>  [33] shape_1.4.6                 DEoptimR_1.0-10            
#>  [35] SparseM_1.81                scales_1.1.1               
#>  [37] DBI_1.1.2                   Rcpp_1.0.8.2               
#>  [39] viridisLite_0.4.0           glmnet_4.1-3               
#>  [41] htmlwidgets_1.5.4           httr_1.4.2                 
#>  [43] ellipsis_0.3.2              pkgconfig_2.0.3            
#>  [45] XML_3.99-0.9                iterpc_0.4.2               
#>  [47] dbplyr_2.1.1                utf8_1.2.2                 
#>  [49] crul_1.2.0                  tidyselect_1.1.2           
#>  [51] rlang_1.0.2                 munsell_0.5.0              
#>  [53] cellranger_1.1.0            tools_4.1.2                
#>  [55] cli_3.1.1                   generics_0.1.2             
#>  [57] broom_0.7.12                mathjaxr_1.6-0             
#>  [59] evaluate_0.15               fastmap_1.1.0              
#>  [61] yaml_2.3.5                  fs_1.5.2                   
#>  [63] zip_2.2.0                   arrangements_1.1.9         
#>  [65] robustbase_0.93-9           nlme_3.1-153               
#>  [67] quantreg_5.88               xml2_1.3.3                 
#>  [69] compiler_4.1.2              rstudioapi_0.13            
#>  [71] plotly_4.10.0               curl_4.3.2                 
#>  [73] reprex_2.0.1                stringi_1.7.6              
#>  [75] gdtools_0.2.4               lattice_0.20-45            
#>  [77] Matrix_1.4-2                nloptr_2.0.0               
#>  [79] vctrs_0.3.8                 CompQuadForm_1.4.3         
#>  [81] pillar_1.7.0                lifecycle_1.0.1            
#>  [83] data.table_1.14.2           bitops_1.0-7               
#>  [85] R6_2.5.1                    BiocIO_1.4.0               
#>  [87] codetools_0.2-18            boot_1.3-28                
#>  [89] MASS_7.3-54                 assertthat_0.2.1           
#>  [91] SummarizedExperiment_1.24.0 rjson_0.2.21               
#>  [93] withr_2.5.0                 nortest_1.0-4              
#>  [95] httpcode_0.3.0              GenomicAlignments_1.30.0   
#>  [97] metafor_3.0-2               Rsamtools_2.10.0           
#>  [99] GenomeInfoDbData_1.2.7      parallel_4.1.2             
#> [101] hms_1.1.1                   grid_4.1.2                 
#> [103] minqa_1.2.4                 rmarkdown_2.13             
#> [105] snakecase_0.11.0            MatrixGenerics_1.6.0       
#> [107] Biobase_2.54.0              lubridate_1.8.0            
#> [109] base64enc_0.1-3             restfulr_0.0.13
```
